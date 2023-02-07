#!/usr/bin/env python3
"""
Genome Analysis Toolkit

https://gatk.broadinstitute.org/
"""
import logging
import os
from typing import Iterable, List, Union

import pysam

import paleomix.common.rtools as rtools
from paleomix.common.command import (
    AtomicCmd,
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


class ApplyBQSRNode(CommandNode):
    def __init__(
        self,
        in_node,
        out_bam,
        out_log=None,
        options={},
        java_options=(),
        dependencies=(),
    ):
        if not isinstance(in_node, BaseRecalibratorNode):
            raise ValueError(in_node)

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
            description="recalibrating base qualities for {!r}".format(in_node.in_bam),
            dependencies=dependencies,
        )


class ApplyVQSRNode(CommandNode):
    def __init__(
        self,
        mode,
        in_vcf,
        in_node,
        out_vcf,
        out_log=None,
        options={},
        java_options=(),
        dependencies=(),
    ):
        if not isinstance(in_node, VariantRecalibratorNode):
            raise ValueError(in_node)
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
            description="recalibrating {}s for {!r}".format(mode, in_vcf),
            dependencies=dependencies,
        )


class BaseRecalibratorNode(CommandNode):
    def __init__(
        self,
        in_reference,
        in_known_sites,
        in_bam,
        out_table,
        out_log=None,
        options={},
        java_options=(),
        dependencies=(),
    ):
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
            description="training base recalibrator on {!r}".format(in_bam),
            dependencies=dependencies,
        )


class CombineGVCFsNode(CommandNode):
    def __init__(
        self,
        in_reference,
        in_variants,
        out_vcf,
        out_log=None,
        options={},
        java_options=(),
        dependencies=(),
    ):
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
            description="combining {} VCFs".format(len(in_variants)),
            dependencies=dependencies,
        )


class CreateSequenceDictionaryNode(CommandNode):
    def __init__(self, in_fasta, options={}, java_options=(), dependencies=()):
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
            description="creating sequence dictionary for {!r}".format(in_fasta),
            dependencies=dependencies,
        )


class FastqToSamNode(CommandNode):
    def __init__(
        self,
        in_fastq,
        out_bam,
        options={},
        java_options=(),
        dependencies=(),
    ) -> None:
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
            description="converting {} to BAM".format(in_fastq),
            dependencies=dependencies,
        )

    def _setup(self, temp: PathTypes) -> None:
        self._enabled = False
        with open_rb(self._in_fastq) as handle:
            # Empty files cannot be processed by GATK
            self._enabled = bool(handle.readline().strip())

        # Prerequisites should always be checked, even if command is not run
        return CommandNode._setup(self, temp)

    def _run(self, temp: PathTypes) -> None:
        if self._enabled:
            return CommandNode._run(self, temp)

        return Node._run(self, temp)

    def _teardown(self, temp: PathTypes) -> None:
        if self._enabled:
            return CommandNode._teardown(self, temp)

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
        in_vcfs,
        out_vcf,
        options={},
        java_options=(),
        dependencies=(),
    ) -> None:
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
            description="gathering {} VCFs into {!r}".format(len(in_vcfs), out_vcf),
            dependencies=dependencies,
        )


class GenotypeGVCFs(CommandNode):
    def __init__(
        self,
        in_reference,
        in_gvcf,
        out_vcf,
        out_log=None,
        options={},
        java_options=(),
        dependencies=(),
    ):
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
            description="genotyping {!r}".format(in_gvcf),
            dependencies=dependencies,
        )


class HaplotypeCallerNode(CommandNode):
    def __init__(
        self,
        in_reference,
        in_bam,
        out_vcf,
        options={},
        out_log=None,
        java_options=(),
        dependencies=(),
    ):
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
            description="calling haplotypes for {!r}".format(in_bam),
            threads=options.get("--native-pair-hmm-threads", 4),
            dependencies=dependencies,
        )


class SplitIntervalsNode(CommandNode):
    def __init__(
        self,
        in_reference: str,
        out_folder: str,
        scatter_count: int = 1,
        options: OptionsType = {},
        java_options: Iterable[str] = (),
        dependencies: Iterable[Node] = (),
    ):
        if scatter_count < 1:
            raise ValueError("scatter_count must be >= 1, not {}".format(scatter_count))

        self.intervals = []
        for idx in range(scatter_count):
            name = "{:03}_of_{:03}".format(idx, scatter_count)
            filename = os.path.join(out_folder, "{:04}.interval_list".format(idx))

            # A list of dicts is used so that sort order is always predictable; this
            # is needed for some operations that require intervals in genome order
            self.intervals.append({"name": name, "filename": filename})

        extra_files: List[Union[InputFile, OutputFile]] = []
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
            description="splitting {!r} into {} intervals".format(
                in_reference, scatter_count
            ),
            dependencies=dependencies,
        )


class TranchesPlotsNode(CommandNode):
    """Task for plotting tranche statistics generated by VariantRecalibrator using a
    custom R-script. This is needed due to the built-in script producing poor output
    in some cases, and to make plotting independant of the variant calling task.
    """

    def __init__(self, input_table, output_prefix, dependencies=()):
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
            description="plotting recalibration tranches for {!r}".format(input_table),
            dependencies=dependencies,
        )


class VariantRecalibratorNode(CommandNode):
    def __init__(
        self,
        mode,
        in_reference,
        in_variant,
        out_recal,
        out_tranches,
        out_r_plot,
        out_log=None,
        options={},
        java_options={},
        dependencies=(),
    ):
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
            description="training {} recalibrator for {!r}".format(mode, in_variant),
            dependencies=dependencies,
        )


def _gatk_command(
    tool,
    java_options,
    tool_options,
    user_tool_options,
    stdout=None,
    stderr=None,
    extra_files=(),
):
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


def _normalize_idx_extension(command, out_bam):
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
