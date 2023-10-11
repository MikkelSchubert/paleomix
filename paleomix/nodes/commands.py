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
"""Implements nodes for calling PALEOMIX commands.

Each node is equivalent to a particular command:
    $ paleomix [...]
"""
from __future__ import annotations

import os
from typing import Iterable

from paleomix.common.command import InputFile, OptionsType, OutputFile, ParallelCmds
from paleomix.common.fileutils import describe_files
from paleomix.node import CommandNode, Node
from paleomix.nodes.samtools import merge_bam_files_command
from paleomix.tools import factory


class CoverageNode(CommandNode):
    def __init__(
        self,
        *,
        target_name: str,
        input_file: str,
        output_file: str,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = factory.new(
            [
                "coverage",
                InputFile(input_file),
                OutputFile(output_file),
                "--target-name",
                target_name,
            ]
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"calculating coverage for {input_file}",
            dependencies=dependencies,
        )


class DepthHistogramNode(CommandNode):
    def __init__(
        self,
        *,
        target_name: str,
        input_file: str,
        output_file: str,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = factory.new(
            [
                "depths",
                InputFile(input_file),
                OutputFile(output_file),
                "--target-name",
                target_name,
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"calculating depth histogram for {input_file}",
            dependencies=dependencies,
        )


class FilterCollapsedBAMNode(CommandNode):
    def __init__(
        self,
        *,
        input_bams: Iterable[str],
        output_bam: str,
        keep_dupes: bool = True,
        dependencies: Iterable[Node] = (),
    ) -> None:
        input_bams = tuple(input_bams)
        if len(input_bams) > 1:
            merge = merge_bam_files_command(input_bams)
            markdup = factory.new(
                "rmdup_collapsed",
                stdin=merge,
                stdout=output_bam,
            )
            command = ParallelCmds([merge, markdup])
        elif len(input_bams) == 1:
            markdup = command = factory.new(
                ["rmdup_collapsed", InputFile(input_bams[0])],
                stdout=output_bam,
            )
        else:
            raise ValueError("empty list of BAM files")

        if not keep_dupes:
            markdup.append("--remove-duplicates")

        CommandNode.__init__(
            self,
            command=command,
            description="filtering merged-read PCR duplicates in {}".format(
                describe_files(input_bams)
            ),
            dependencies=dependencies,
        )


class FinalizeBAMNode(CommandNode):
    def __init__(
        self,
        in_bams: Iterable[str],
        out_passed: str,
        out_failed: str,
        out_json: str,
        threads: int = 1,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        in_bams = tuple(in_bams)
        if len(in_bams) > 1:
            merge = merge_bam_files_command(in_bams)
            finalize = factory.new(
                "ngs:finalize_bam",
                stdin=merge,
            )

            command = ParallelCmds([merge, finalize])
        elif len(in_bams) == 1:
            finalize = command = factory.new(
                ["ngs:finalize_bam", InputFile(in_bams[0])]
            )
        else:
            raise ValueError(in_bams)

        finalize.merge_options(
            user_options=options,
            fixed_options={
                "--out-passed": OutputFile(out_passed),
                "--out-failed": OutputFile(out_failed),
                "--out-json": OutputFile(out_json),
                "--threads": threads,
            },
        )

        CommandNode.__init__(
            self,
            description="creating finalized BAMs {!r} and {!r}".format(
                os.path.basename(out_passed),
                os.path.basename(out_failed),
            ),
            command=command,
            dependencies=dependencies,
            threads=threads,
        )
