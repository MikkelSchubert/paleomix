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

from paleomix.node import CommandNode
from paleomix.atomiccmd.builder import (
    AtomicCmdBuilder,
    apply_options,
)

import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions


_VERSION_CHECK = versions.Requirement(
    call=("AdapterRemoval", "--version"),
    search=r"ver. (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(2, 2, 0),
)


class SE_AdapterRemovalNode(CommandNode):
    def __init__(
        self, input_file, output_prefix, threads=1, options={}, dependencies=()
    ):
        # See below for parameters in common between SE/PE
        cmd = _get_common_parameters(threads=threads, options=options)

        # Prefix for output files, ensure that all end up in temp folder
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        output_tmpl = output_prefix + ".%s.gz"
        cmd.set_kwargs(
            TEMP_OUT_BASENAME=os.path.basename(output_prefix),
            OUT_SETTINGS=output_prefix + ".settings",
            OUT_MATE_1=output_tmpl % ("truncated",),
            OUT_DISCARDED=output_tmpl % ("discarded",),
        )

        cmd.set_option("--file1", "%(IN_READS_1)s")
        cmd.set_kwargs(IN_READS_1=input_file)

        apply_options(cmd, options)

        CommandNode.__init__(
            self,
            command=cmd.finalize(),
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
        cmd = _get_common_parameters(threads=threads, options=options)

        # Prefix for output files, to ensure that all end up in temp folder
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        output_tmpl = output_prefix + ".%s.gz"
        cmd.set_kwargs(
            TEMP_OUT_BASENAME=os.path.basename(output_prefix),
            OUT_SETTINGS=output_prefix + ".settings",
            OUT_READS_1=output_tmpl % ("pair1.truncated",),
            OUT_READS_2=output_tmpl % ("pair2.truncated",),
            OUT_SINGLETON=output_tmpl % ("singleton.truncated",),
            OUT_DISCARDED=output_tmpl % ("discarded",),
        )

        if collapse:
            cmd.set_option("--collapse")

            cmd.set_kwargs(
                OUT_COLLAPSED=output_tmpl % ("collapsed",),
                OUT_COLLAPSED_TRUNC=output_tmpl % ("collapsed.truncated",),
            )

        cmd.set_option("--file1", "%(IN_READS_1)s")
        cmd.set_option("--file2", "%(IN_READS_2)s")
        cmd.set_kwargs(IN_READS_1=input_file_1, IN_READS_2=input_file_2)

        apply_options(cmd, options)

        CommandNode.__init__(
            self,
            command=cmd.finalize(),
            threads=threads,
            description="trimming PE adapters from %s"
            % fileutils.describe_paired_files(input_file_1, input_file_2),
            dependencies=dependencies,
        )


def _get_common_parameters(options, threads=1):
    cmd = AtomicCmdBuilder("AdapterRemoval", CHECK_VERSION=_VERSION_CHECK)

    # Gzip compress FASTQ files
    cmd.set_option("--gzip")

    # Trim Ns at read ends
    cmd.set_option("--trimns", fixed=False)
    # Trim low quality scores
    cmd.set_option("--trimqualities", fixed=False)

    # Fix number of threads to ensure consistency when scheduling node
    cmd.set_option("--threads", threads)

    # Ensure that any user-specified list of adapters is tracked
    adapter_list = options.pop("--adapter-list", None)
    if adapter_list is not None:
        cmd.cmd.set_option("--adapter-list", "%(IN_ADAPTER_LIST)s")
        cmd.command.set_kwargs(IN_ADAPTER_LIST=adapter_list)

    for key in ("--trim5p", "--trim3p"):
        values = options.pop(key, None)
        if values is not None:
            cmd.add_value(key)
            if isinstance(values, list):
                for value in values:
                    cmd.add_value(value)
            else:
                cmd.add_value(values)

    return cmd
