#!/usr/bin/env python3
"""
FastQC - A quality control analysis tool for high throughput sequencing data

https://github.com/s-andrews/FastQC
"""
import os
import re

from paleomix.atomiccmd.command import AtomicCmd, InputFile, OutputFile
from paleomix.common.versions import Requirement
from paleomix.node import CommandNode


# File extensions striped by FASTQ for output filenames
_FASTQC_EXCLUDED_EXTENSIONS = re.compile(
    r"(\.gz|\.bz2|\.txt|\.fastq|\.fq|\.csfastq|\.sam|\.bam)+$"
)


class FastQCNode(CommandNode):
    def __init__(self, in_file, out_folder, options={}, dependencies=()):
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
                    search=r"FastQC v(\d+).(\d+).(\d+)",
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
            description="fastQC of {}".format(in_file),
            dependencies=dependencies,
        )
