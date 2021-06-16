#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
from typing import Any, Iterable, Type, Union

import paleomix.common.versions as versions
from paleomix.common.command import (
    AtomicCmd,
    InputFile,
    OptionsType,
    OutputFile,
    ParallelCmds,
    TempOutputFile,
)
from paleomix.node import CommandNode, Node, NodeError
from paleomix.nodes.bwa import (
    _get_max_threads,
    _get_node_description,
    _new_cleanup_command,
)

BOWTIE2_VERSION = versions.Requirement(
    call=("bowtie2", "--version"),
    search=r"version (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(2, 3, 0),
)


class Bowtie2IndexNode(CommandNode):
    def __init__(self, input_file: str, dependencies: Iterable[Node] = ()):
        command = _bowtie2_template(
            (
                "bowtie2-build",
                InputFile(input_file),
                TempOutputFile(input_file),
            ),
            reference=input_file,
            iotype=OutputFile,
        )

        CommandNode.__init__(
            self,
            command=command,
            description="creating Bowtie2 index for %s" % (input_file,),
            dependencies=dependencies,
        )


class Bowtie2Node(CommandNode):
    def __init__(
        self,
        input_file_1: str,
        input_file_2: str,
        output_file: str,
        reference: str,
        threads: int = 2,
        mapping_options: OptionsType = {},
        cleanup_options: OptionsType = {},
        dependencies: Iterable[Node] = (),
    ):
        aln = _bowtie2_template(
            ["bowtie2"],
            reference=reference,
            stdout=AtomicCmd.PIPE,
        )

        fixed_options: OptionsType = {
            "--threads": _get_max_threads(reference, threads),
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

        cleanup = _new_cleanup_command(
            stdin=aln,
            in_reference=reference,
            out_bam=output_file,
            max_threads=fixed_options["--threads"],
            paired_end=input_file_1 and input_file_2,
            options=cleanup_options,
        )

        description = _get_node_description(
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


def _bowtie2_template(
    call: Any,
    reference: str,
    iotype: Union[Type[InputFile], Type[OutputFile]] = InputFile,
    **kwargs: Any
):
    return AtomicCmd(
        call,
        extra_files=[
            iotype(reference + postfix)
            for postfix in (
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            )
        ],
        requirements=[BOWTIE2_VERSION],
        **kwargs
    )
