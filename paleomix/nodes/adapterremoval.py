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
import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions

from paleomix.node import CommandNode
from paleomix.atomiccmd.command import AtomicCmd, InputFile, OutputFile, TempOutputFile


_VERSION_CHECK = versions.Requirement(
    call=("AdapterRemoval", "--version"),
    search=r"ver. (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(2, 2, 0),
)


class SE_AdapterRemovalNode(CommandNode):
    def __init__(
        self, input_file, output_prefix, threads=1, options={}, dependencies=()
    ):
        self.out_settings = output_prefix + ".settings"
        self.out_truncated = "{}.truncated.gz".format(output_prefix)
        self.out_discarded = "{}.discarded.gz".format(output_prefix)

        command = AtomicCmd(
            "AdapterRemoval",
            extra_files=[
                OutputFile(self.out_settings),
                OutputFile(self.out_truncated),
                OutputFile(self.out_discarded),
            ],
            requirements=[_VERSION_CHECK],
        )

        # Ignored for SE reads
        options = dict(options)
        options.pop("--collapse", None)
        options.pop("--collapse-deterministic", None)
        options.pop("--collapse-conservatively", None)

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
                "--basename": TempOutputFile(output_prefix),
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
        threads=1,
        options={},
        dependencies=(),
    ):
        self.out_settings = output_prefix + ".settings"
        self.out_paired = "{}.pair{{Pair}}.truncated.gz".format(output_prefix)
        self.out_singleton = "{}.singleton.truncated.gz".format(output_prefix)
        self.out_discarded = "{}.discarded.gz".format(output_prefix)
        self.out_merged = None
        self.out_merged_truncated = None

        command = AtomicCmd(
            "AdapterRemoval",
            extra_files=[
                OutputFile(self.out_settings),
                OutputFile(self.out_paired.format(Pair=1)),
                OutputFile(self.out_paired.format(Pair=2)),
                OutputFile(self.out_singleton),
                OutputFile(self.out_discarded),
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
            "--basename": TempOutputFile(output_prefix),
        }

        if options.keys() & (
            "--collapse",
            "--collapse-deterministic",
            "--collapse-conservatively",
        ):
            self.out_merged = "{}.collapsed.gz".format(output_prefix)
            self.out_merged_truncated = "{}.collapsed.truncated.gz".format(
                output_prefix
            )

            command.add_extra_files(
                [
                    OutputFile(self.out_merged),
                    OutputFile(self.out_merged_truncated),
                ]
            )

        # Ensure that any user-specified list of adapters is tracked
        if "--adapter-list" in options:
            options = dict(options)
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
