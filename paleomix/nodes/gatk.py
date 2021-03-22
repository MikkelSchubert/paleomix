#!/usr/bin/env python3
"""
Genome Analysis Toolkit

https://gatk.broadinstitute.org/
"""
import os

import paleomix.common.system

from paleomix.atomiccmd.sets import SequentialCmds
from paleomix.atomiccmd.command2 import AtomicCmd2, InputFile, OutputFile
from paleomix.common.fileutils import swap_ext
from paleomix.node import CommandNode, NodeError


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

        self.out_bam = out_bam

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
                "--variant": InputFile(in_node.in_variant),
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
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="recalibrating {}s for {!r}".format(mode, in_node.in_variant),
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
        self.out_vcf = out_vcf

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
    ):
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


class GatherVcfsNode(CommandNode):
    def __init__(
        self,
        in_vcfs,
        out_vcf,
        options={},
        java_options=(),
        dependencies=(),
    ):
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


class MarkDuplicatesNode(CommandNode):
    def __init__(
        self,
        in_bams,
        out_bam,
        out_metrics=None,
        options={},
        java_options=(),
        dependencies=(),
    ):
        out_metrics = out_metrics or swap_ext(out_bam, ".metrics")

        command = _gatk_command(
            tool="MarkDuplicates",
            tool_options={
                "--OUTPUT": OutputFile(out_bam),
                "--METRICS_FILE": OutputFile(out_metrics),
                # FIXME: Workaround for CSI indexed BAM files
                # Validation is mostly left to manual ValidateSamFile runs; required
                # because .csi indexed BAM records can have "invalid" bins.
                "--VALIDATION_STRINGENCY": "LENIENT",
            },
            user_tool_options=options,
            java_options=java_options,
        )

        _set_max_open_files(command, "--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP")

        for filename in in_bams:
            command.append("--INPUT", InputFile(filename))

        CommandNode.__init__(
            self,
            command=command,
            description="detecting PCR duplicates in {} BAM files".format(len(in_bams)),
            dependencies=dependencies,
        )


class MergeSamFilesNode(CommandNode):
    def __init__(self, in_bams, out_bam, options={}, java_options=(), dependencies=()):
        command = _gatk_command(
            tool="MergeSamFiles",
            tool_options={
                "--INPUT": [InputFile(filename) for filename in in_bams],
                "--OUTPUT": OutputFile(out_bam),
                "--SORT_ORDER": "coordinate",
                # FIXME: Workaround for CSI indexed BAM files
                # Validation is mostly left to manual ValidateSamFile runs; required
                # because .csi indexed BAM records can have "invalid" bins.
                "--VALIDATION_STRINGENCY": "LENIENT",
            },
            user_tool_options=options,
            java_options=java_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            description="merging {} file(s) into {}".format(len(in_bams), out_bam),
            dependencies=dependencies,
        )


class ValidateBAMNode(CommandNode):
    def __init__(
        self,
        in_bam,
        in_index=None,
        out_log=None,
        big_genome_mode=False,
        options={},
        java_options=(),
        dependencies=(),
    ):
        self.out_log = out_log or swap_ext(in_bam, ".validated")
        # Allows the node to be used as if it was the mapping/indexing node itself
        self.out_bam = in_bam
        self.out_idx = in_index

        command = _gatk_command(
            tool="ValidateSamFile",
            tool_options={
                "--INPUT": InputFile(in_bam),
                "--OUTPUT": OutputFile(self.out_log),
            },
            user_tool_options=options,
            java_options=java_options,
        )

        # FIXME: Workaround for check with high rate of false positives
        # Ignored due to rate of false positives for runs with few hits,
        # where high-quality reads may cause mis-identification of qualities
        command.append("--IGNORE", "INVALID_QUALITY_FORMAT")

        _set_max_open_files(command, "--MAX_OPEN_TEMP_FILES")

        if big_genome_mode:
            self._configure_for_big_genome(command)

        if in_index:
            command.add_extra_files([InputFile(in_index)])

        CommandNode.__init__(
            self,
            command=command,
            description="validating {!r}".format(in_bam),
            dependencies=dependencies,
        )

    def _configure_for_big_genome(self, command):
        # FIXME: Workaround for CSI indexed BAM files
        # Validation is mostly left to manual ValidateSamFile runs; required
        # because .csi indexed BAM records can have "invalid" bins.
        command.append("--IGNORE", "INVALID_INDEXING_BIN")

        # FIXME: Workaround for useless warning; no BAI indexes for large genomes
        command.append("--IGNORE", "REF_SEQ_TOO_LONG_FOR_BAI")


class SplitIntervalsNode(CommandNode):
    def __init__(
        self,
        in_reference,
        out_folder,
        scatter_count=1,
        options={},
        java_options=(),
        dependencies=(),
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

        extra_files = [OutputFile(interval["filename"]) for interval in self.intervals]
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


class VariantRecalibratorNode(CommandNode):
    def __init__(
        self,
        mode,
        in_reference,
        in_variant,
        out_prefix,
        options={},
        out_log=None,
        java_options={},
        dependencies=(),
    ):
        if mode not in ("INDEL", "SNP"):
            raise ValueError(mode)

        self.in_reference = in_reference
        self.in_variant = in_variant
        self.out_recal = out_prefix + ".recal.vcf.gz"
        self.out_tranches = out_prefix + ".tranches"
        self.out_r_plot = out_prefix + ".r"
        self.out_log = out_prefix + ".log"

        options = dict(options)
        extra_files = [
            InputFile(in_reference + ".fai"),
            InputFile(swap_ext(in_reference, ".dict")),
            # FIXME: Extension depends on fname: .idx for .vcf, .tbi for vcf.gz, ..
            OutputFile(self.out_recal + ".tbi"),
            OutputFile(self.out_r_plot + ".pdf"),
        ]

        if mode == "SNP":
            # This plot is not generated in INDEL mode
            extra_files.append(OutputFile(self.out_tranches + ".pdf"))

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
            stderr=out_log,
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
    command = AtomicCmd2(
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

    command.add_extra_files([OutputFile(in_index, temporary=True)])

    rename = AtomicCmd2(
        [
            "mv",
            InputFile(in_index, temporary=True),
            OutputFile(out_bam + ".bai"),
        ]
    )

    return SequentialCmds([command, rename])


# Maximum number of open files
_MAX_OPEN_FILES = None
# Fraction of per-process max open files to use
_FRAC_MAX_OPEN_FILES = 0.95
# Default maximum number of open temporary files used by Picard
_DEFAULT_MAX_OPEN_FILES = 8000


def _set_max_open_files(command, key):
    """Sets the maximum number of open files a picard process should use, at most.
    Conservatively lowered than the actual ulimit.
    """
    global _MAX_OPEN_FILES
    if _MAX_OPEN_FILES is None:
        _MAX_OPEN_FILES = (paleomix.common.system.get_max_open_files(),)

    (max_open_files,) = _MAX_OPEN_FILES
    if max_open_files:
        max_open_files = int(max_open_files * _FRAC_MAX_OPEN_FILES)

        if max_open_files < _DEFAULT_MAX_OPEN_FILES:
            command.append(key, max_open_files)
