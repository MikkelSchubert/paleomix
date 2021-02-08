#!/usr/bin/env python3
"""
Genome Analysis Toolkit

https://gatk.broadinstitute.org/
"""
import os

from paleomix.atomiccmd.sets import SequentialCmds
from paleomix.atomiccmd.command2 import AtomicCmd2, InputFile, OutputFile, CmdError
from paleomix.common.fileutils import swap_ext
from paleomix.common.utilities import fill_dict
from paleomix.node import CommandNode


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
        )

        # Rename ${FILENAME}.bai to ${FILENAME}.bam.bai
        command = _normalize_idx_extension(command, out_bam)

        CommandNode.__init__(
            self,
            command=command,
            description="recalibrating base qualities",
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
                # FIXME: Index depends on vcf extension
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="recalibrating %ss for %r" % (mode, in_node.in_variant),
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
                InputFile(swap_ext(in_reference, ".dict")),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="training base recalibrator",
            dependencies=dependencies,
        )


class CombineGVCFs(CommandNode):
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
                # FIXME: Index depends on vcf extension
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="genotyping %r" % (in_gvcf,),
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
                # FIXME: Index depends on vcf extension
                OutputFile(out_vcf + ".tbi"),
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="calling haplotypes",
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
        resource_files = []
        for key, value in options.items():
            if key.startswith("--resource"):
                options[key] = InputFile(value)
                # FIXME: Expected index may also depend on extension (see above)
                resource_files.append(InputFile(value + ".tbi"))

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
            extra_files=[
                # FIXME: Extension depends on fname: .idx for .vcf, .tbi for vcf.gz, ..
                OutputFile(self.out_recal + ".tbi"),
                OutputFile(self.out_tranches + ".pdf"),
                OutputFile(self.out_r_plot + ".pdf"),
            ]
            + resource_files,
        )

        CommandNode.__init__(
            self,
            command=command,
            description="training %s recalibrator for %r" % (mode, in_variant),
            dependencies=dependencies,
        )


def _gatk_command(
    tool,
    java_options,
    tool_options,
    user_tool_options,
    stderr=None,
    extra_files=(),
):
    command = AtomicCmd2(["gatk"], stderr=stderr, extra_files=extra_files)
    if java_options:
        command.append("--java-options", " ".join(java_options))
    command.append(tool)

    options = dict(tool_options)
    for key, value in user_tool_options.items():
        if key in options:
            raise CmdError("cannot set %r for GATK %s" % (key, tool))

        options[key] = value

    command.append_options(options)

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
