#!/usr/bin/env python3
"""
fastp - An ultra-fast all-in-one FASTQ preprocessor

https://github.com/OpenGene/fastp
"""
from paleomix.atomiccmd.command2 import AtomicCmd2, InputFile, OutputFile
from paleomix.common.fileutils import describe_paired_files, reroot_path
from paleomix.node import CommandNode, NodeError


_NOT_SUPPORTED = set(["--overlapped_out", "--stdin", "--stdout", "--interleaved_in"])


class FastpNode(CommandNode):
    def __init__(
        self,
        in_fq_1,
        in_fq_2,
        out_fq_1,
        out_fq_2,
        out_merged,
        out_unpaired,
        out_html,
        out_json,
        options={},
        dependencies=(),
    ):
        for key in _NOT_SUPPORTED & set(dict(options)):
            raise NotImplementedError("fastp option {} not supported".format(key))

        command = AtomicCmd2(
            ["fastp"],
        )

        # All unpaired reads are written to the same file
        out_unpaired = OutputFile(out_unpaired)

        fastp_options = {
            "--merge": None,
            "--in1": InputFile(in_fq_1),
            "--in2": InputFile(in_fq_2),
            "--out1": OutputFile(out_fq_1),
            "--out2": OutputFile(out_fq_2),
            "--merged_out": OutputFile(out_merged),
            "--unpaired1": out_unpaired,
            "--unpaired2": out_unpaired,
            "--html": OutputFile(out_html),
            "--json": OutputFile(out_json),
        }

        for key in set(options) & set(fastp_options):
            raise NodeError("cannot set override fastp option {}".format(key))

        fastp_options.update(options)
        command.append_options(fastp_options)

        CommandNode.__init__(
            self,
            command=command,
            description="pre-processing %s"
            % (describe_paired_files([in_fq_1], [in_fq_2]),),
            threads=options.get("--thread", 2),
            dependencies=dependencies,
        )
