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
MultiQC - Aggregate results from bioinformatics analyses across many samples into a
single report

https://multiqc.info/
"""
from __future__ import annotations

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
