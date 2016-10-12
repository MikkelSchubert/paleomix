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
    CommandNode, \
    NodeError
from paleomix.atomiccmd.command import \
    AtomicCmd
from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder, \
    use_customizable_cli_parameters, \
    create_customizable_cli_parameters
from paleomix.atomiccmd.sets import \
    ParallelCmds
from paleomix.nodes.bwa import \
    _get_node_description, \
    _process_output, \
    _get_max_threads
from paleomix.common.utilities import \
    safe_coerce_to_tuple

import paleomix.common.versions as versions


BOWTIE2_VERSION = versions.Requirement(call=("bowtie2", "--version"),
                                       search=r"version (\d+)\.(\d+)\.(\d+)",
                                       checks=versions.GE(2, 1, 0))


class Bowtie2IndexNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file, prefix=None, dependencies=()):
        prefix = prefix if prefix else input_file
        params = _bowtie2_template(("bowtie2-build"), prefix, iotype="OUT",
                                   IN_FILE=input_file,
                                   TEMP_OUT_PREFIX=os.path.basename(prefix),
                                   CHECK_VERSION=BOWTIE2_VERSION)
        params.add_value("%(IN_FILE)s")
        # Destination prefix, in temp folder
        params.add_value("%(TEMP_OUT_PREFIX)s")

        return {"prefix": prefix,
                "command": params,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()
        description = "<Bowtie2 Index '%s' -> '%s.*'>" \
            % (parameters.input_file, parameters.prefix)

        CommandNode.__init__(self,
                             command=command,
                             description=description,
                             dependencies=parameters.dependencies)


class Bowtie2Node(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file_1, input_file_2, output_file, reference,
                  prefix, threads=2, log_file=None, dependencies=()):

        # Setting IN_FILE_2 to None makes AtomicCmd ignore this key
        aln = _bowtie2_template(("bowtie2",), prefix,
                                OUT_STDOUT=AtomicCmd.PIPE,
                                CHECK_VERSION=BOWTIE2_VERSION)
        aln.set_option("-x", prefix)

        if log_file is not None:
            aln.set_kwargs(OUT_STDERR=log_file)

        if input_file_1 and not input_file_2:
            aln.add_multiple_options("-U", safe_coerce_to_tuple(input_file_1),
                                     template="IN_FILE_1_%02i")
        elif input_file_1 and input_file_2:
            aln.add_multiple_options("-1", safe_coerce_to_tuple(input_file_1),
                                     template="IN_FILE_1_%02i")
            aln.add_multiple_options("-2", safe_coerce_to_tuple(input_file_2),
                                     template="IN_FILE_2_%02i")
        else:
            raise NodeError("Input 1, OR both input 1 and input 2 must "
                            "be specified for Bowtie2 node")

        max_threads = _get_max_threads(reference, threads)
        aln.set_option("--threads", max_threads)

        run_fixmate = input_file_1 and input_file_2
        order, commands = _process_output(aln, output_file, reference,
                                          run_fixmate=run_fixmate)
        commands["aln"] = aln

        return {"commands": commands,
                "order": ["aln"] + order,
                "threads": max_threads,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = ParallelCmds([parameters.commands[key].finalize()
                               for key in parameters.order])

        algorithm = "PE" if parameters.input_file_2 else "SE"
        description \
            = _get_node_description(name="Bowtie2",
                                    algorithm=algorithm,
                                    input_files_1=parameters.input_file_1,
                                    input_files_2=parameters.input_file_2,
                                    prefix=parameters.prefix,
                                    threads=parameters.threads)

        CommandNode.__init__(self,
                             command=command,
                             description=description,
                             threads=parameters.threads,
                             dependencies=parameters.dependencies)


def _bowtie2_template(call, prefix, iotype="IN", **kwargs):
    for postfix in ("1.bt2", "2.bt2", "3.bt2", "4.bt2",
                    "rev.1.bt2", "rev.2.bt2"):
        key = "%s_PREFIX_%s" % (iotype, postfix.upper())
        kwargs[key] = prefix + "." + postfix

    return AtomicCmdBuilder(call, **kwargs)
