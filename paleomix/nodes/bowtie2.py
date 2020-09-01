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

from paleomix.node import CommandNode, NodeError
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.builder import (
    AtomicCmdBuilder,
    apply_options,
)
from paleomix.atomiccmd.sets import ParallelCmds
from paleomix.nodes.bwa import (
    _get_node_description,
    _new_cleanup_command,
    _get_max_threads,
)

import paleomix.common.versions as versions


BOWTIE2_VERSION = versions.Requirement(
    call=("bowtie2", "--version"),
    search=r"version (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(2, 3, 0),
)


class Bowtie2IndexNode(CommandNode):
    def __init__(self, input_file, prefix=None, dependencies=()):
        prefix = prefix if prefix else input_file
        builder = _bowtie2_template(
            ("bowtie2-build"),
            prefix,
            iotype="OUT",
            IN_FILE=input_file,
            TEMP_OUT_PREFIX=os.path.basename(prefix),
            CHECK_VERSION=BOWTIE2_VERSION,
        )

        builder.add_value("%(IN_FILE)s")
        # Destination prefix, in temp folder
        builder.add_value("%(TEMP_OUT_PREFIX)s")

        CommandNode.__init__(
            self,
            command=builder.finalize(),
            description="creating Bowtie2 index for %s" % (input_file,),
            dependencies=dependencies,
        )


class Bowtie2Node(CommandNode):
    def __init__(
        self,
        input_file_1,
        input_file_2,
        output_file,
        reference,
        prefix,
        threads=2,
        log_file=None,
        mapping_options={},
        cleanup_options={},
        dependencies=(),
    ):
        # Setting IN_FILE_2 to None makes AtomicCmd ignore this key
        aln = _bowtie2_template(
            ("bowtie2",),
            prefix,
            OUT_STDOUT=AtomicCmd.PIPE,
            CHECK_VERSION=BOWTIE2_VERSION,
        )

        aln.set_option("-x", prefix)

        if log_file is not None:
            aln.set_kwargs(OUT_STDERR=log_file)

        if input_file_1 and not input_file_2:
            aln.add_option("-U", input_file_1)
        elif input_file_1 and input_file_2:
            aln.add_option("-1", input_file_1)
            aln.add_option("-2", input_file_2)
        else:
            raise NodeError(
                "Input 1, OR both input 1 and input 2 must "
                "be specified for Bowtie2 node"
            )

        max_threads = _get_max_threads(reference, threads)
        aln.set_option("--threads", max_threads)

        cleanup = _new_cleanup_command(
            aln, output_file, reference, paired_end=input_file_1 and input_file_2
        )

        apply_options(aln, mapping_options)
        apply_options(cleanup, cleanup_options)

        algorithm = "PE" if input_file_2 else "SE"
        description = _get_node_description(
            name="Bowtie2",
            algorithm=algorithm,
            input_files_1=input_file_1,
            input_files_2=input_file_2,
            prefix=prefix,
            threads=threads,
        )

        CommandNode.__init__(
            self,
            command=ParallelCmds([aln.finalize(), cleanup.finalize()]),
            description=description,
            threads=threads,
            dependencies=dependencies,
        )


def _bowtie2_template(call, prefix, iotype="IN", **kwargs):
    for postfix in ("1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"):
        key = "%s_PREFIX_%s" % (iotype, postfix.upper())
        kwargs[key] = prefix + "." + postfix

    return AtomicCmdBuilder(call, **kwargs)
