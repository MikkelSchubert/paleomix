#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
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

from paleomix.node import \
    CommandNode
from paleomix.atomiccmd.sets import \
    ParallelCmds
from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder, \
    use_customizable_cli_parameters, \
    create_customizable_cli_parameters

import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions
import paleomix.tools.factory as factory


_VERSION_CHECK = versions.Requirement(call=("AdapterRemoval", "--version"),
                                      search=r"ver. (\d+)\.(\d+)\.(\d+)",
                                      checks=versions.GE(2, 1, 5))


class SE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_files, output_prefix, output_format="bz2",
                  threads=1, dependencies=()):
        # See below for parameters in common between SE/PE
        cmd = _get_common_parameters(output_format, threads=threads)

        # Prefix for output files, ensure that all end up in temp folder
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        output_tmpl = output_prefix + ".%s." + output_format
        cmd.set_kwargs(TEMP_OUT_BASENAME=os.path.basename(output_prefix),

                       OUT_SETTINGS=output_prefix + ".settings",
                       OUT_MATE_1=output_tmpl % ("truncated",),
                       OUT_DISCARDED=output_tmpl % ("discarded",))

        if len(input_files) > 1:
            # Uncompressed reads (piped from 'paleomix cat')
            cmd.set_option("--file1", "%(TEMP_IN_READS_1)s")
            cmd.set_kwargs(TEMP_IN_READS_1="uncompressed_input")
        else:
            cmd.set_option("--file1", "%(IN_READS_1)s")
            cmd.set_kwargs(IN_READS_1=input_files[0])

        return {"command": cmd,
                "threads": threads,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()

        self._multi_file_input = len(parameters.input_files) > 1
        if self._multi_file_input:
            cat = _build_cat_command(parameters.input_files, "uncompressed_input")
            command = ParallelCmds((command, cat))

        CommandNode.__init__(self,
                             command=command,
                             threads=parameters.threads,
                             description="<AdapterRM (SE): %s -> '%s.*'>"
                             % (fileutils.describe_files(parameters.input_files),
                                parameters.output_prefix),
                             dependencies=parameters.dependencies)

    def _setup(self, config, temp):
        if self._multi_file_input:
            os.mkfifo(os.path.join(os.path.join(temp, "uncompressed_input")))

        CommandNode._setup(self, config, temp)


class PE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_files_1, input_files_2, output_prefix,
                  output_format="bz2", collapse=True, threads=1,
                  dependencies=()):
        if len(input_files_1) != len(input_files_2):
            raise ValueError("Unequal number of mate 1 and mate 2 files")

        cmd = _get_common_parameters(output_format, threads=threads)

        # Prefix for output files, to ensure that all end up in temp folder
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        output_tmpl = output_prefix + ".%s." + output_format
        cmd.set_kwargs(TEMP_OUT_BASENAME=os.path.basename(output_prefix),

                       OUT_SETTINGS=output_prefix + ".settings",
                       OUT_READS_1=output_tmpl % ("pair1.truncated",),
                       OUT_READS_2=output_tmpl % ("pair2.truncated",),
                       OUT_SINGLETON=output_tmpl % ("singleton.truncated",),
                       OUT_DISCARDED=output_tmpl % ("discarded",))

        if collapse:
            cmd.set_option("--collapse")

            cmd.set_kwargs(OUT_COLLAPSED=output_tmpl % ("collapsed",),
                           OUT_COLLAPSED_TRUNC=output_tmpl
                           % ("collapsed.truncated",))

        if len(input_files_1) > 1:
            # Uncompressed reads (piped from 'paleomix cat')
            cmd.set_option("--file1", "%(TEMP_IN_READS_1)s")
            cmd.set_option("--file2", "%(TEMP_IN_READS_2)s")
            cmd.set_kwargs(TEMP_IN_READS_1="uncompressed_input_1",
                           TEMP_IN_READS_2="uncompressed_input_2")
        else:
            cmd.set_option("--file1", "%(IN_READS_1)s")
            cmd.set_option("--file2", "%(IN_READS_2)s")
            cmd.set_kwargs(IN_READS_1=input_files_1[0],
                           IN_READS_2=input_files_2[0])

        return {"command": cmd,
                "threads": threads,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()
        self._multi_file_input = len(parameters.input_files_1) > 1
        if self._multi_file_input:
            cat_1 = _build_cat_command(parameters.input_files_1, "uncompressed_input_1")
            cat_2 = _build_cat_command(parameters.input_files_2, "uncompressed_input_2")
            command = ParallelCmds((command, cat_1, cat_2))

        CommandNode.__init__(self,
                             command=command,
                             threads=parameters.threads,
                             description="<AdapterRM (PE): %s -> '%s.*'>"
                             % (fileutils.describe_paired_files(parameters.input_files_1,
                                                                parameters.input_files_2),
                                parameters.output_prefix),
                             dependencies=parameters.dependencies)

    def _setup(self, config, temp):
        if self._multi_file_input:
            os.mkfifo(os.path.join(os.path.join(temp, "uncompressed_input_1")))
            os.mkfifo(os.path.join(os.path.join(temp, "uncompressed_input_2")))

        CommandNode._setup(self, config, temp)


def _build_cat_command(input_files, output_file):
    cat = factory.new("cat")
    cat.set_option("--output", "%(TEMP_OUT_CAT)s")
    cat.set_kwargs(TEMP_OUT_CAT=output_file)
    cat.add_multiple_values(input_files)

    return cat.finalize()


def _get_common_parameters(output_format, threads=1):
    cmd = AtomicCmdBuilder("AdapterRemoval",
                           CHECK_VERSION=_VERSION_CHECK)

    if output_format == "bz2":
        cmd.set_option("--bzip2")
    elif output_format == "gz":
        cmd.set_option("--gzip")
    else:
        raise ValueError("Invalid output compression %r" % (output_format,))

    # Trim Ns at read ends
    cmd.set_option("--trimns", fixed=False)
    # Trim low quality scores
    cmd.set_option("--trimqualities", fixed=False)

    # Fix number of threads to ensure consistency when scheduling node
    cmd.set_option("--threads", threads)

    return cmd
