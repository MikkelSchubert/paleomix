# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
"""
MultiQC - Aggregate results from bioinformatics analyses across many samples into a
single report

https://multiqc.info/
"""

from __future__ import annotations

import fnmatch
import os
from collections.abc import Iterable
from typing import Literal

from paleomix.common.command import AtomicCmd, InputFile, OptionsType, OutputFile
from paleomix.node import CommandNode, Node, NodeError

# Supported modules and their expected files
# For more, see https://multiqc.info/docs/#multiqc-modules
MODULES = {
    "fastp": "*fastp.json",
    "fastqc": "*_fastqc.zip",
}


class MultiQCNode(CommandNode):
    def __init__(
        self,
        source: Literal["fastp", "fastqc"],
        output_prefix: str,
        dependencies: Iterable[Node],
        options: OptionsType | None = None,
    ) -> None:
        pattern = MODULES.get(source)
        if pattern is None:
            raise NodeError(f"unsupported MultiQC source {source!r}")

        if options is None:
            options = {}

        command = AtomicCmd(
            [
                "multiqc",
                "--zip-data-dir",
                "--filename",
                os.path.basename(output_prefix),
                "--outdir",
                "%(TEMP_DIR)s",
                "--module",
                source,
            ],
            extra_files=(
                OutputFile(output_prefix + ".html"),
                OutputFile(output_prefix + "_data.zip"),
            ),
        )

        command.append_options(options)

        input_files_found = False
        dependencies = list(dependencies)
        for node in dependencies:
            for filename in fnmatch.filter(node.output_files, pattern):
                input_files_found = True
                command.append(InputFile(filename))

        if not input_files_found:
            raise NodeError(f"no {source!r} input files found for MultiQC")

        CommandNode.__init__(
            self,
            command=command,
            description=f"multiQC of {len(dependencies)} reports",
            dependencies=dependencies,
        )
