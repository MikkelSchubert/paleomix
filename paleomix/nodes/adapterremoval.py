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
from paleomix.atomiccmd.builder import (
    AtomicCmdBuilder,
    apply_options,
)
from paleomix.node import CommandNode

_VERSION_2_CHECK = versions.Requirement(
    call=("AdapterRemoval", "--version"),
    search=r"ver. (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(2, 2, 0),
)

_VERSION_3_CHECK = versions.Requirement(
    call=("adapterremoval3", "--version"),
    search=r"v(\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(3, 0, 0),
)


class SE_AdapterRemoval2Node(CommandNode):
    def __init__(
        self, input_file, output_prefix, threads=1, options={}, dependencies=()
    ):
        # See below for parameters in common between SE/PE
        cmd = AtomicCmdBuilder("AdapterRemoval", CHECK_VERSION=_VERSION_2_CHECK)
        cmd = _set_common_parameters(cmd, threads=threads, options=options)

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

        self.files = {
            "Statistics": output_prefix + ".settings",
            "Single": output_tmpl % ("truncated",),
        }

        CommandNode.__init__(
            self,
            command=cmd.finalize(),
            threads=threads,
            description="trimming SE adapters from %s"
            % fileutils.describe_files(input_file),
            dependencies=dependencies,
        )


class PE_AdapterRemoval2Node(CommandNode):
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
        # See below for parameters in common between SE/PE
        cmd = AtomicCmdBuilder("AdapterRemoval", CHECK_VERSION=_VERSION_2_CHECK)
        cmd = _set_common_parameters(cmd, threads=threads, options=options)

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

        self.files = {
            "Statistics": output_prefix + ".settings",
            "Singleton": output_tmpl % ("singleton.truncated",),
            "Paired": output_tmpl % ("pair{Pair}.truncated",),
        }

        if collapse:
            cmd.set_option("--collapse")
            cmd.set_kwargs(
                OUT_COLLAPSED=output_tmpl % ("collapsed",),
                OUT_COLLAPSED_TRUNC=output_tmpl % ("collapsed.truncated",),
            )

            self.files["Collapsed"] = output_tmpl % ("collapsed",)
            self.files["CollapsedTruncated"] = output_tmpl % ("collapsed.truncated",)

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


class SE_AdapterRemoval3Node(CommandNode):
    def __init__(
        self,
        input_file,
        output_prefix,
        threads=1,
        options=None,
        dependencies=(),
    ):
        options = {} if options is None else dict(options)
        # Option is no longer supported (output is always Phred+33)
        options.pop("--qualitybase-output", None)

        # See below for parameters in common between SE/PE
        cmd = AtomicCmdBuilder("adapterremoval3", CHECK_VERSION=_VERSION_3_CHECK)
        cmd = _set_common_parameters(cmd, threads=threads, options=options)

        # Prefix for output files, ensure that all end up in temp folder
        cmd.set_option("--out-prefix", "%(TEMP_OUT_BASENAME)s")

        output_tmpl = output_prefix + ".%s.fastq.gz"
        cmd.set_kwargs(
            TEMP_OUT_BASENAME=os.path.basename(output_prefix),
            OUT_JSON=output_prefix + ".json",
            OUT_HTML=output_prefix + ".html",
            OUT_MATE_1=output_tmpl % ("r1",),
            OUT_DISCARDED=output_tmpl % ("discarded",),
            IN_READS_1=input_file,
        )

        cmd.set_option("--in-file1", "%(IN_READS_1)s")
        cmd.set_option("--out-discarded", "%(OUT_DISCARDED)s")

        apply_options(cmd, options)

        self.files = {
            "Statistics": output_prefix + ".json",
            "Single": output_tmpl % ("r1",),
        }

        CommandNode.__init__(
            self,
            command=cmd.finalize(),
            threads=threads,
            description="trimming SE adapters from %s"
            % fileutils.describe_files(input_file),
            dependencies=dependencies,
        )


class PE_AdapterRemoval3Node(CommandNode):
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
        options = {} if options is None else dict(options)
        # Option is no longer supported (output is always Phred+33)
        options.pop("--qualitybase-output", None)

        # See below for parameters in common between SE/PE
        cmd = AtomicCmdBuilder("adapterremoval3", CHECK_VERSION=_VERSION_3_CHECK)
        cmd = _set_common_parameters(cmd, threads=threads, options=options)

        # Prefix for output files, to ensure that all end up in temp folder
        cmd.set_option("--out-prefix", "%(TEMP_OUT_BASENAME)s")

        output_tmpl = output_prefix + ".%s.fastq.gz"
        cmd.set_kwargs(
            TEMP_OUT_BASENAME=os.path.basename(output_prefix),
            OUT_JSON=output_prefix + ".json",
            OUT_HTML=output_prefix + ".html",
            OUT_READS_1=output_tmpl % ("r1",),
            OUT_READS_2=output_tmpl % ("r2",),
            OUT_SINGLETON=output_tmpl % ("singleton",),
            OUT_DISCARDED=output_tmpl % ("discarded",),
        )

        self.files = {
            "Statistics": output_prefix + ".json",
            "Singleton": output_tmpl % ("singleton",),
            "Paired": output_tmpl % ("r{Pair}",),
        }

        if collapse:
            cmd.set_option("--merge")
            cmd.set_kwargs(OUT_MERGED=output_tmpl % ("merged",))
            self.files["Collapsed"] = output_tmpl % ("merged",)

        cmd.set_option("--in-file1", "%(IN_READS_1)s")
        cmd.set_option("--in-file2", "%(IN_READS_2)s")
        cmd.set_kwargs(IN_READS_1=input_file_1, IN_READS_2=input_file_2)
        cmd.set_option("--out-discarded", "%(OUT_DISCARDED)s")

        apply_options(cmd, options)

        CommandNode.__init__(
            self,
            command=cmd.finalize(),
            threads=threads,
            description="trimming PE adapters from %s"
            % fileutils.describe_paired_files(input_file_1, input_file_2),
            dependencies=dependencies,
        )


def _set_common_parameters(cmd, options, threads=1):
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
        cmd.set_option("--adapter-list", "%(IN_ADAPTER_LIST)s")
        cmd.set_kwargs(IN_ADAPTER_LIST=adapter_list)

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
