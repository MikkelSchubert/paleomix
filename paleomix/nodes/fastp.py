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
"""
fastp - An ultra-fast all-in-one FASTQ preprocessor

https://github.com/OpenGene/fastp
"""
from __future__ import annotations

from typing import Iterable

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
