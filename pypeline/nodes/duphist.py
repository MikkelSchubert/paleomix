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
from pypeline.atomiccmd.command import \
    AtomicCmd
from pypeline.atomiccmd.sets import \
    ParallelCmds
from pypeline.common.fileutils import \
    describe_files
from pypeline.nodes.picard import \
    MultiBAMInput, \
    MultiBAMInputNode


class DuplicateHistogramNode(MultiBAMInputNode):
    def __init__(self, config, input_files, output_file, dependencies=()):
        bam_input = MultiBAMInput(config, input_files)
        duphist_command = AtomicCmd(['bam_duphist', '%(TEMP_IN_BAM)s'],
                                    TEMP_IN_BAM=bam_input.pipe,
                                    OUT_STDOUT=output_file)
        commands = ParallelCmds(bam_input.commands + [duphist_command])

        description = "<DuplicateHistogram: %s -> %r>" \
            % (describe_files(input_files), output_file)
        MultiBAMInputNode.__init__(self,
                                   bam_input=bam_input,
                                   command=commands,
                                   description=description,
                                   dependencies=dependencies)
