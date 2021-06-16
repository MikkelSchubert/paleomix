#!/usr/bin/env python3
"""
MultiQC - Aggregate results from bioinformatics analyses across many samples into a
single report

https://multiqc.info/
"""
import fnmatch
import os

from paleomix.common.command import AtomicCmd, InputFile, OutputFile
from paleomix.node import CommandNode, NodeError

# Supported modules and their expected files
# For more, see https://multiqc.info/docs/#multiqc-modules
MODULES = {
    "fastp": "*fastp.json",
    "fastqc": "*_fastqc.zip",
}


class MultiQCNode(CommandNode):
    def __init__(self, source, output_prefix, dependencies, options={}):
        pattern = MODULES.get(source)
        if pattern is None:
            raise NodeError("unsupported MultiQC source {!r}".format(source))

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
        for node in dependencies:
            for filename in fnmatch.filter(node.output_files, pattern):
                input_files_found = True
                command.append(InputFile(filename))

        if not input_files_found:
            raise NodeError("no {!r} input files found for MultiQC".format(source))

        CommandNode.__init__(
            self,
            command=command,
            description="multiQC of %i reports" % (len(dependencies),),
            dependencies=dependencies,
        )
