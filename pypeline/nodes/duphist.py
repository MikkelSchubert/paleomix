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
import random
import collections

import pysam

from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
from pypeline.common.fileutils import \
    describe_files, \
    move_file
from pypeline.nodes.picard import \
    concatenate_input_bams


class DuplicateHistogramNode(CommandNode):
    def __init__(self, config, input_files, output_file, dependencies=()):
        cat_cmds, cat_obj = concatenate_input_bams(config, input_files)
        duphist_command = AtomicCmd(['bam_duphist', '-'],
                                    IN_STDIN=cat_obj,
                                    OUT_STDOUT=output_file)
        commands = ParallelCmds(cat_cmds + [duphist_command])

        description = "<DuplicateHistogram: %s -> %r>" \
            % (describe_files(input_files), output_file)
        CommandNode.__init__(self,
                             command=commands,
                             description=description,
                             dependencies=dependencies)
