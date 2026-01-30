# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

from collections.abc import Iterable
from typing import Literal

from paleomix.common.command import InputFile
from paleomix.common.fileutils import describe_files
from paleomix.node import CommandNode, Node
from paleomix.tools import factory


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
