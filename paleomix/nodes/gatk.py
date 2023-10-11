#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING, Iterable

import pysam

from paleomix.common import rtools
from paleomix.common.command import (
    AtomicCmd,
    AtomicFileTypes,
    AuxiliaryFile,
    InputFile,
    OptionsType,
    OutputFile,
    SequentialCmds,
    TempInputFile,
    TempOutputFile,
)
from paleomix.common.fileutils import (
    PathTypes,
    move_file,
    open_rb,
    reroot_path,
    swap_ext,
)
from paleomix.node import CommandNode, Node, NodeError

if TYPE_CHECKING:
    from typing_extensions import Literal


class ApplyBQSRNode(CommandNode):
    def __init__(
        self,
        in_node: BaseRecalibratorNode,
        out_bam: str,
        out_log: str | None = None,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        if not isinstance(in_node, BaseRecalibratorNode):
            raise TypeError(in_node)

        # WORKAROUND: ApplyBQSR defaults to using compression level 2 for output,
        #             resulting in significantly larger BAMs. Revert that behavior
        #             unless the user has already done so. GATK 4.2.3.0 and more.
        if not any(o.startswith("-Dsamjdk.compression_level=") for o in java_options):
            java_options = list(java_options)
            java_options.append("-Dsamjdk.compression_level=6")

        command = _gatk_command(
            tool="ApplyBQSR",
            tool_options={
                "--reference": InputFile(in_node.in_reference),
                "--input": InputFile(in_node.in_bam),
                "--bqsr-recal-file": InputFile(in_node.out_table),
                "--output": OutputFile(out_bam),
            },
            user_tool_options=options,
            java_options=java_options,
            stderr=out_log,
            extra_files=[
                InputFile(in_node.in_bam + ".bai"),
                InputFile(in_node.in_reference + ".fai"),
                InputFile(swap_ext(in_node.in_reference, ".dict")),
            ],
        )

        # Rename ${FILENAME}.bai to ${FILENAME}.bam.bai
        command = _normalize_idx_extension(command, out_bam)

        CommandNode.__init__(
            self,
            command=command,
            description=f"recalibrating base qualities for {in_node.in_bam!r}",
            dependencies=dependencies,
        )


class ApplyVQSRNode(CommandNode):
    def __init__(
        self,
        mode: Literal["INDEL", "SNP"],
        in_vcf: str,
        in_node: str,
        out_vcf: str,
        out_log: str | None = None,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        if not isinstance(in_node, VariantRecalibratorNode):
            raise TypeError(in_node)
        elif mode not in ("INDEL", "SNP"):
            raise ValueError(mode)

        self.out_vcf = out_vcf

        command = _gatk_command(
            tool="ApplyVQSR",
            tool_options={
                "--mode": mode,
                "--reference": InputFile(in_node.in_reference),
                "--variant": InputFile(in_vcf),
                "--recal-file": InputFile(in_node.out_recal),
                "--tranches-file": InputFile(in_node.out_tranches),
                "--output": OutputFile(out_vcf),
            },
            user_tool_options=options,
            java_options=java_options,
            stderr=out_log,
            extra_files=[
                InputFile(in_node.in_reference + ".fai"),
                InputFile(swap_ext(in_node.in_reference, ".dict")),
                # FIXME: Index depends on vcf extension
                InputFile(in_vcf + ".tbi"),
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"recalibrating {mode}s for {in_vcf!r}",
            dependencies=dependencies,
        )


class BaseRecalibratorNode(CommandNode):
    def __init__(
        self,
        in_reference: str,
        in_known_sites: str,
        in_bam: str,
        out_table: str,
        out_log: str | None = None,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        self.in_reference = in_reference
        self.in_bam = in_bam
        self.out_table = out_table

        command = _gatk_command(
            tool="BaseRecalibrator",
            tool_options={
                "--reference": InputFile(in_reference),
                "--known-sites": InputFile(in_known_sites),
                "--input": InputFile(in_bam),
                "--output": OutputFile(out_table),
            },
            user_tool_options=options,
            java_options=java_options,
            stderr=out_log,
            extra_files=[
                InputFile(in_bam + ".bai"),
                InputFile(in_reference + ".fai"),
                InputFile(swap_ext(in_reference, ".dict")),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"training base recalibrator on {in_bam!r}",
            dependencies=dependencies,
        )


class CombineGVCFsNode(CommandNode):
    def __init__(
        self,
        in_reference: str,
        in_variants: str,
        out_vcf: str,
        out_log: str | None = None,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        command = _gatk_command(
            tool="CombineGVCFs",
            tool_options={
                "--reference": InputFile(in_reference),
                "--variant": [InputFile(fname) for fname in in_variants],
                "--output": OutputFile(out_vcf),
            },
            user_tool_options=options,
            java_options=java_options,
            stderr=out_log,
            extra_files=[
                InputFile(in_reference + ".fai"),
                InputFile(swap_ext(in_reference, ".dict")),
                # FIXME: Index depends on vcf extension
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        # FIXME: Are TBI files actually required? Index also depends on format
        command.add_extra_files(InputFile(fname + ".tbi") for fname in in_variants)

        CommandNode.__init__(
            self,
            command=command,
            description=f"combining {len(in_variants)} VCFs",
            dependencies=dependencies,
        )


class CreateSequenceDictionaryNode(CommandNode):
    def __init__(
        self,
        *,
        in_fasta: str,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        command = _gatk_command(
            tool="CreateSequenceDictionary",
            tool_options={
                "--REFERENCE": InputFile(in_fasta),
                "--OUTPUT": OutputFile(swap_ext(in_fasta, ".dict")),
            },
            user_tool_options=options,
            java_options=java_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"creating sequence dictionary for {in_fasta!r}",
            dependencies=dependencies,
        )


class FastqToSamNode(CommandNode):
    def __init__(
        self,
        in_fastq: str,
        out_bam: str,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        self._enabled = True
        self._in_fastq = in_fastq
        self._out_bam = out_bam
        self._sort_order = options.get("--SORT-ORDER")

        if "--SAMPLE_NAME" not in dict(options):
            raise NodeError("--SAMPLE_NAME must be specified for FastqToSamNode")

        command = _gatk_command(
            tool="FastqToSam",
            tool_options={
                "--FASTQ": InputFile(in_fastq),
                "--OUTPUT": OutputFile(out_bam),
            },
            user_tool_options=options,
            java_options=java_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"converting {in_fastq} to BAM",
            dependencies=dependencies,
        )

    def _setup(self, temp: PathTypes) -> None:
        self._enabled = False
        with open_rb(self._in_fastq) as handle:
            # Empty files cannot be processed by GATK
            self._enabled = bool(handle.readline().strip())

        # Prerequisites should always be checked, even if command is not run
        return super()._setup(temp)

    def _run(self, temp: PathTypes) -> None:
        if self._enabled:
            super()._run(temp)

    def _teardown(self, temp: PathTypes) -> None:
        if self._enabled:
            return super()._teardown(temp)

        log = logging.getLogger(__name__)
        log.debug("creating dummy BAM file for %r", self._out_bam)

        header = {"HD": {"VN": "1.6"}}
        if self._sort_order:
            header["HD"]["SO"] = self._sort_order

        filename = reroot_path(temp, self._out_bam)
        with pysam.AlignmentFile(filename, "wb", header=header):
            pass

        move_file(filename, self._out_bam)

        return Node._teardown(self, temp)


class GatherVcfsNode(CommandNode):
    def __init__(
        self,
        in_vcfs: str,
        out_vcf: str,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        command = _gatk_command(
            tool="GatherVcfs",
            tool_options={
                "--INPUT": [InputFile(in_vcf) for in_vcf in in_vcfs],
                "--OUTPUT": OutputFile(out_vcf),
            },
            user_tool_options=options,
            java_options=java_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"gathering {len(in_vcfs)} VCFs into {out_vcf!r}",
            dependencies=dependencies,
        )


class GenotypeGVCFs(CommandNode):
    def __init__(
        self,
        *,
        in_reference: str,
        in_gvcf: str,
        out_vcf: str,
        out_log: str | None = None,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        self.out_vcf = out_vcf

        command = _gatk_command(
            tool="GenotypeGVCFs",
            tool_options={
                "--reference": InputFile(in_reference),
                "--variant": InputFile(in_gvcf),
                "--output": OutputFile(out_vcf),
            },
            user_tool_options=options,
            java_options=java_options,
            stderr=out_log,
            extra_files=[
                InputFile(in_reference + ".fai"),
                InputFile(swap_ext(in_reference, ".dict")),
                # FIXME: Index depends on vcf extension
                InputFile(in_gvcf + ".tbi"),
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"genotyping {in_gvcf!r}",
            dependencies=dependencies,
        )


class HaplotypeCallerNode(CommandNode):
    def __init__(
        self,
        *,
        in_reference: str,
        in_bam: str,
        out_vcf: str,
        out_log: str | None = None,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        self.out_vcf = out_vcf

        command = _gatk_command(
            tool="HaplotypeCaller",
            tool_options={
                "--reference": InputFile(in_reference),
                "--input": InputFile(in_bam),
                "--output": OutputFile(out_vcf),
            },
            user_tool_options=options,
            java_options=java_options,
            stderr=out_log,
            extra_files=[
                InputFile(in_reference + ".fai"),
                InputFile(swap_ext(in_reference, ".dict")),
                # FIXME: Index depends on vcf extension
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"calling haplotypes for {in_bam!r}",
            threads=options.get("--native-pair-hmm-threads", 4),
            dependencies=dependencies,
        )


class SplitIntervalsNode(CommandNode):
    def __init__(
        self,
        in_reference: str,
        out_folder: str,
        scatter_count: int = 1,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        if scatter_count < 1:
            raise ValueError(f"scatter_count must be >= 1, not {scatter_count}")

        self.intervals: list[dict[str, str]] = []
        for idx in range(scatter_count):
            name = f"{idx:03}_of_{scatter_count:03}"
            filename = os.path.join(out_folder, f"{idx:04}.interval_list")

            # A list of dicts is used so that sort order is always predictable; this
            # is needed for some operations that require intervals in genome order
            self.intervals.append({"name": name, "filename": filename})

        extra_files: list[InputFile | OutputFile] = []
        for interval in self.intervals:
            extra_files.append(OutputFile(interval["filename"]))

        # Both the fai and dict files are required by the command
        extra_files.append(InputFile(in_reference + ".fai"))
        extra_files.append(InputFile(swap_ext(in_reference, ".dict")))

        command = _gatk_command(
            tool="SplitIntervals",
            tool_options={
                "--reference": InputFile(in_reference),
                "--output": "%(TEMP_DIR)s",
                "--scatter-count": scatter_count,
                "--extension": ".interval_list",
            },
            user_tool_options=options,
            java_options=java_options,
            extra_files=extra_files,
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"splitting {in_reference!r} into {scatter_count} intervals",
            dependencies=dependencies,
        )


class TranchesPlotsNode(CommandNode):
    """Task for plotting tranche statistics generated by VariantRecalibrator using a
    custom R-script. This is needed due to the built-in script producing poor output
    in some cases, and to make plotting independant of the variant calling task.
    """

    def __init__(
        self,
        input_table: str,
        output_prefix: str,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = AtomicCmd(
            (
                "Rscript",
                AuxiliaryFile(rtools.rscript("ngs", "tranches.r")),
                InputFile(input_table),
                TempOutputFile(output_prefix),
            ),
            extra_files=(
                OutputFile(output_prefix + ".pdf"),
                OutputFile(output_prefix + ".txt"),
            ),
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"plotting recalibration tranches for {input_table!r}",
            dependencies=dependencies,
        )


class VariantRecalibratorNode(CommandNode):
    def __init__(
        self,
        *,
        mode: Literal["INDEL", "SNP"],
        in_reference: str,
        in_variant: str,
        out_recal: str,
        out_tranches: str,
        out_r_plot: str,
        out_log: str | None = None,
        options: OptionsType | None = None,
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        if mode not in ("INDEL", "SNP"):
            raise ValueError(mode)

        self.in_reference = in_reference
        self.in_variant = in_variant
        self.out_recal = out_recal
        self.out_tranches = out_tranches
        self.out_r_plot = out_r_plot
        self.out_log = out_log

        options = dict(options)
        extra_files = [
            InputFile(in_reference + ".fai"),
            InputFile(swap_ext(in_reference, ".dict")),
            # FIXME: Extension depends on fname: .idx for .vcf, .tbi for vcf.gz, ..
            InputFile(self.in_variant + ".tbi"),
            OutputFile(self.out_recal + ".tbi"),
            OutputFile(self.out_r_plot + ".pdf"),
        ]

        # WORKAROUND: GATK creates unreadable tranche plots due to unsorted data, so we
        #             use a custom R-script instead (via TranchesPlotsNode). This issue
        #             applies to GATK 4.2.4.1 and others.
        extra_files.append(TempOutputFile(self.out_tranches + ".pdf"))

        for key, value in options.items():
            if key.startswith("--resource"):
                if not isinstance(value, str):
                    raise TypeError(f"--resource={value!r} is not a string")

                options[key] = InputFile(value)
                # FIXME: Expected index may also depend on extension (see above)
                extra_files.append(InputFile(value + ".tbi"))

        command = _gatk_command(
            tool="VariantRecalibrator",
            tool_options={
                "--mode": mode,
                "--reference": InputFile(in_reference),
                "--variant": InputFile(in_variant),
                "--output": OutputFile(self.out_recal),
                "--tranches-file": OutputFile(self.out_tranches),
                "--rscript-file": OutputFile(self.out_r_plot),
            },
            user_tool_options=options,
            java_options=java_options,
            stderr=self.out_log,
            extra_files=extra_files,
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"training {mode} recalibrator for {in_variant!r}",
            dependencies=dependencies,
        )


def _gatk_command(
    tool: str,
    java_options: Iterable[str],
    tool_options: OptionsType,
    user_tool_options: OptionsType,
    stdout: int | str | OutputFile | None = None,
    stderr: int | str | OutputFile | None = None,
    extra_files: Iterable[AtomicFileTypes] = (),
) -> AtomicCmd:
    command = AtomicCmd(
        ["gatk"],
        stdout=stdout,
        stderr=stderr,
        extra_files=extra_files,
    )

    if java_options:
        command.append("--java-options", " ".join(java_options))

    command.append(tool)
    command.merge_options(
        user_options=user_tool_options,
        fixed_options=tool_options,
    )

    return command


def _normalize_idx_extension(command: AtomicCmd, out_bam: str) -> SequentialCmds:
    # FIXME: Support for other index formats?
    in_index = swap_ext(os.path.basename(out_bam), ".bai")

    command.add_extra_files([TempOutputFile(in_index)])

    rename = AtomicCmd(
        [
            "mv",
            TempInputFile(in_index),
            OutputFile(out_bam + ".bai"),
        ]
    )

    return SequentialCmds([command, rename])
