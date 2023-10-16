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

from typing import TYPE_CHECKING, Iterable

from paleomix.common.command import InputFile
from paleomix.common.fileutils import describe_files
from paleomix.node import CommandNode, Node
from paleomix.tools import factory

if TYPE_CHECKING:
    from typing_extensions import Literal


class ValidateFASTQFilesNode(CommandNode):
    def __init__(
        self,
        *,
        input_files: Iterable[str],
        output_file: str,
        offset: Literal[33, 64, "solexa"],
        collapsed: bool = False,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = factory.new(
            [":validate_fastq", "--offset", str(offset)],
            stdout=output_file,
        )

        if collapsed:
            command.append("--collapsed")

        for filename in input_files:
            command.append(InputFile(filename))

        CommandNode.__init__(
            self,
            description=f"validating {describe_files(input_files)}",
            command=command,
            dependencies=dependencies,
        )


class ValidateFASTAFilesNode(CommandNode):
    def __init__(
        self,
        input_file: str,
        output_file: str,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = factory.new(
            [":validate_fasta", InputFile(input_file)],
            stdout=output_file,
        )

        CommandNode.__init__(
            self,
            description=f"validating {input_file}",
            command=command,
            dependencies=dependencies,
        )
