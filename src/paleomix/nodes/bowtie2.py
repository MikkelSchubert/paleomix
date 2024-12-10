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

from collections.abc import Iterable

from paleomix.common import versions
from paleomix.common.command import (
    AtomicCmd,
    AtomicFileTypes,
    InputFile,
    OptionsType,
    OutputFile,
    ParallelCmds,
    TempOutputFile,
)
from paleomix.node import CommandNode, Node, NodeError
from paleomix.nodes.bwa import (
    get_max_threads,
    get_node_description,
    new_cleanup_command,
)

BOWTIE2_VERSION = versions.Requirement(
    call=("bowtie2", "--version"),
    regexp=r"version (\d+\.\d+\.\d+)",
    specifiers=">=2.3.0",
)


class Bowtie2IndexNode(CommandNode):
    def __init__(self, *, input_file: str, dependencies: Iterable[Node] = ()) -> None:
        command = AtomicCmd(
            (
                "bowtie2-build",
                InputFile(input_file),
                TempOutputFile(input_file),
            ),
            extra_files=_reference_files(input_file, iotype=OutputFile),
            requirements=[BOWTIE2_VERSION],
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"creating Bowtie2 index for {input_file}",
            dependencies=dependencies,
        )


class Bowtie2Node(CommandNode):
    def __init__(
        self,
        *,
        input_file_1: str,
        input_file_2: str | None,
        output_file: str,
        reference: str,
        threads: int = 2,
        mapping_options: OptionsType | None = None,
        cleanup_options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if cleanup_options is None:
            cleanup_options = {}
        if mapping_options is None:
            mapping_options = {}

        aln = AtomicCmd(
            ["bowtie2"],
            extra_files=_reference_files(reference, iotype=InputFile),
            stdout=AtomicCmd.PIPE,
            requirements=[BOWTIE2_VERSION],
        )

        threads = get_max_threads(reference, threads)
        fixed_options: OptionsType = {
            "--threads": threads,
            "-x": reference,
        }

        if input_file_1 and not input_file_2:
            fixed_options["-U"] = input_file_1
        elif input_file_1 and input_file_2:
            fixed_options["-1"] = input_file_1
            fixed_options["-2"] = input_file_2
        else:
            raise NodeError(
                "Input 1, OR both input 1 and input 2 must "
                "be specified for Bowtie2 node"
            )

        aln.merge_options(
            user_options=mapping_options,
            fixed_options=fixed_options,
        )

        cleanup = new_cleanup_command(
            stdin=aln,
            in_reference=reference,
            out_bam=output_file,
            max_threads=threads,
            paired_end=bool(input_file_1 and input_file_2),
            options=cleanup_options,
        )

        description = get_node_description(
            name="Bowtie2",
            input_files_1=input_file_1,
            input_files_2=input_file_2,
            reference=reference,
        )

        CommandNode.__init__(
            self,
            command=ParallelCmds([aln, cleanup]),
            description=description,
            threads=threads,
            dependencies=dependencies,
        )


def _reference_files(
    reference: str,
    iotype: type[InputFile | OutputFile],
) -> Iterable[AtomicFileTypes]:
    for postfix in (
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    ):
        yield iotype(reference + postfix)
