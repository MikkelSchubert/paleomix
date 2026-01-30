# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
"""
FastQC - A quality control analysis tool for high throughput sequencing data

https://github.com/s-andrews/FastQC
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable

from paleomix.common.command import AtomicCmd, InputFile, OptionsType, OutputFile
from paleomix.common.versions import Requirement
from paleomix.node import CommandNode, Node

# File extensions striped by FASTQ for output filenames
_FASTQC_EXCLUDED_EXTENSIONS = re.compile(
    r"(\.gz|\.bz2|\.txt|\.fastq|\.fq|\.csfastq|\.sam|\.bam)+$"
)


class FastQCNode(CommandNode):
    def __init__(
        self,
        in_file: str,
        out_folder: str,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        out_prefix = _FASTQC_EXCLUDED_EXTENSIONS.sub("", os.path.basename(in_file))
        command = AtomicCmd(
            ["fastqc", InputFile(in_file)],
            extra_files=[
                OutputFile(os.path.join(out_folder, out_prefix + "_fastqc.html")),
                OutputFile(os.path.join(out_folder, out_prefix + "_fastqc.zip")),
            ],
            requirements=[
                Requirement(
                    name="FastQC",
                    call=["fastqc", "--version"],
                    regexp=r"FastQC v(\d+\.\d+\.\d+)",
                ),
            ],
        )

        command.merge_options(
            user_options=options,
            fixed_options={
                "--outdir": "%(TEMP_DIR)s",
                "--dir": "%(TEMP_DIR)s",
            },
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"fastQC of {in_file}",
            dependencies=dependencies,
        )
