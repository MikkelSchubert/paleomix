# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
"""
fastp - An ultra-fast all-in-one FASTQ preprocessor

https://github.com/OpenGene/fastp
"""

from __future__ import annotations

from collections.abc import Iterable

from paleomix.common.command import AtomicCmd, InputFile, OptionsType, OutputFile
from paleomix.common.fileutils import describe_paired_files
from paleomix.node import CommandNode, Node


class FastpNode(CommandNode):
    def __init__(
        self,
        *,
        in_fq_1: str,
        in_fq_2: str,
        out_fq_1: str,
        out_fq_2: str,
        out_merged: str,
        out_unpaired: str,
        out_failed: str,
        out_html: str,
        out_json: str,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        command = AtomicCmd(
            ["fastp"],
        )

        # All unpaired reads are written to the same file
        out_unpaired_file = OutputFile(out_unpaired)

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
                "--unpaired1": out_unpaired_file,
                "--unpaired2": out_unpaired_file,
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

        threads = options.get("--thread", 2)
        if not isinstance(threads, int):
            raise TypeError(threads)

        CommandNode.__init__(
            self,
            command=command,
            description=f"pre-processing {describe_paired_files([in_fq_1], [in_fq_2])}",
            # FIXME: Number of threads seems hard to predict; about --thread +1/+2?
            threads=threads + 1,
            dependencies=dependencies,
        )
