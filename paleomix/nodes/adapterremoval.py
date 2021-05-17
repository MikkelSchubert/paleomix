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
import os

import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions

from paleomix.node import CommandNode
from paleomix.atomiccmd.command2 import AtomicCmd2, InputFile, OutputFile


_VERSION_CHECK = versions.Requirement(
    call=("AdapterRemoval", "--version"),
    search=r"ver. (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(2, 2, 0),
)


class SE_AdapterRemovalNode(CommandNode):
    def __init__(
        self, input_file, output_prefix, threads=1, options={}, dependencies=()
    ):
        command = AtomicCmd2(
            "AdapterRemoval",
            extra_files=[
                OutputFile(output_prefix + ".settings"),
                OutputFile("{}.truncated.gz".format(output_prefix)),
                OutputFile("{}.discarded.gz".format(output_prefix)),
            ],
            requirements=[_VERSION_CHECK],
        )

        # Ensure that any user-specified list of adapters is tracked
        if "--adapter-list" in options:
            options["--adapter-list"] = InputFile(options["--adapter-list"])

        command.merge_options(
            user_options=options,
            fixed_options={
                "--file1": InputFile(input_file),
                # Gzip compress FASTQ files
                "--gzip": None,
                # Fix number of threads to ensure consistency when scheduling node
                "--threads": threads,
                # Prefix for output files, ensure that all end up in temp folder
                "--basename": OutputFile(
                    os.path.basename(output_prefix),
                    temporary=True,
                ),
            },
        )

        CommandNode.__init__(
            self,
            command=command,
            threads=threads,
            description="trimming SE adapters from %s"
            % fileutils.describe_files(input_file),
            dependencies=dependencies,
        )


class PE_AdapterRemovalNode(CommandNode):
    def __init__(
        self,
        input_file_1,
        input_file_2,
        output_prefix,
        collapse=True,
        threads=1,
        options={},
        dependencies=(),
    ):

        command = AtomicCmd2(
            "AdapterRemoval",
            extra_files=[
                OutputFile(output_prefix + ".settings"),
                OutputFile("{}.pair1.truncated.gz".format(output_prefix)),
                OutputFile("{}.pair2.truncated.gz".format(output_prefix)),
                OutputFile("{}.singleton.truncated.gz".format(output_prefix)),
                OutputFile("{}.discarded.gz".format(output_prefix)),
            ],
            requirements=[_VERSION_CHECK],
        )

        fixed_options = {
            "--file1": InputFile(input_file_1),
            "--file2": InputFile(input_file_2),
            # Gzip compress FASTQ files
            "--gzip": None,
            # Fix number of threads to ensure consistency when scheduling node
            "--threads": threads,
            # Prefix for output files, ensure that all end up in temp folder
            "--basename": OutputFile(os.path.basename(output_prefix), temporary=True),
        }

        if collapse:
            fixed_options["--collapse"] = None
            command.add_extra_files(
                [
                    OutputFile("{}.collapsed.gz".format(output_prefix)),
                    OutputFile("{}.collapsed.truncated.gz".format(output_prefix)),
                ]
            )

        # Ensure that any user-specified list of adapters is tracked
        if "--adapter-list" in options:
            options["--adapter-list"] = InputFile(options["--adapter-list"])

        command.merge_options(
            user_options=options,
            fixed_options=fixed_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            threads=threads,
            description="trimming PE adapters from %s"
            % fileutils.describe_paired_files(input_file_1, input_file_2),
            dependencies=dependencies,
        )
