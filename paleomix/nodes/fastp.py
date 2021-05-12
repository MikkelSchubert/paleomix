#!/usr/bin/env python3
"""
fastp - An ultra-fast all-in-one FASTQ preprocessor

https://github.com/OpenGene/fastp
"""
from paleomix.atomiccmd.command2 import AtomicCmd2, InputFile, OutputFile
from paleomix.common.fileutils import describe_paired_files
from paleomix.node import CommandNode


class FastpNode(CommandNode):
    def __init__(
        self,
        in_fq_1,
        in_fq_2,
        out_fq_1,
        out_fq_2,
        out_merged,
        out_unpaired,
        out_failed,
        out_html,
        out_json,
        options={},
        dependencies=(),
    ):
        command = AtomicCmd2(
            ["fastp"],
        )

        # All unpaired reads are written to the same file
        out_unpaired = OutputFile(out_unpaired)

        command.merge_options(
            user_options=options,
            fixed_options={
                # FIXME: Make merging optional
                "--merge": None,
                "--in1": InputFile(in_fq_1),
                "--in2": InputFile(in_fq_2),
                "--out1": OutputFile(out_fq_1),
                "--out2": OutputFile(out_fq_2),
                "--failed_out": OutputFile(out_failed),
                "--merged_out": OutputFile(out_merged),
                # Using the same output path for both --unpaired is supported by fastp
                "--unpaired1": out_unpaired,
                "--unpaired2": out_unpaired,
                "--html": OutputFile(out_html),
                "--json": OutputFile(out_json),
            },
            blacklisted_options=[
                "--overlapped_out",
                "--stdin",
                "--stdout",
                "--interleaved_in",
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="pre-processing %s"
            % (describe_paired_files([in_fq_1], [in_fq_2]),),
            # FIXME: Number of threads seems hard to predict; about --thread +1/+2?
            threads=options.get("--thread", 2) + 1,
            dependencies=dependencies,
        )
