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
FastQC - A quality control analysis tool for high throughput sequencing data

https://github.com/s-andrews/FastQC
"""
from __future__ import annotations

import os
import re
from typing import Iterable

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
        if options is None:
            options = {}

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
