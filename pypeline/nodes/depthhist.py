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
from pypeline.node import \
    CommandNode
from pypeline.atomiccmd.builder import \
    AtomicCmdBuilder


class DepthHistogramNode(CommandNode):
    def __init__(self, target_name, input_file, output_file,
                 regions_file=None, dependencies=()):
        call = ['bam_depths', '%(IN_FILE)s', '%(OUT_FILE)s',
                '--target-name', target_name]
        builder = AtomicCmdBuilder(call,
                                   IN_FILE=input_file,
                                   OUT_FILE=output_file)

        if regions_file:
            builder.set_option('--regions-file', '%(IN_REGIONS)s')
            builder.set_kwargs(IN_REGIONS=regions_file)

        command = builder.finalize()

        CommandNode.__init__(self,
                             command=command,
                             description="<DepthHistogram: %s -> '%s'>"
                                         % (input_file, output_file),
                             dependencies=dependencies)
