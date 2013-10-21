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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd

import pypeline.common.versions as versions


BEDTOOLS_VERSION = versions.Requirement(call   = ("bedtools", "--version"),
                                        search = r"bedtools v?(\d+)\.(\d+)\.(\d+)",
                                        checks = versions.GE(2, 15, 0))




class SlopBedNode(CommandNode):
    def __init__(self, infile, outfile, genome, from_start = 0, from_end = 0, strand_relative = False, dependencies = ()):
        if type(from_start) != type(from_end):
            raise ValueError("'from_start' and 'from_end' should be of same type!")

        call = ["bedtools", "slop",
                "-i", "%(IN_FILE)s",
                "-g", "%(IN_GENOME)s",
                "-l", str(from_start),
                "-r", str(from_end)]

        if strand_relative:
            call.append("-s")
        if type(from_start) is float:
            call.append("-pct")

        command = AtomicCmd(call,
                            IN_FILE    = infile,
                            IN_GENOME  = genome,
                            OUT_STDOUT = outfile,
                            CHECK_VERSION = BEDTOOLS_VERSION)

        description = "<SlopBed: '%s' -> '%s'>" % (infile, outfile)

        CommandNode.__init__(self,
                             description  = description,
                             command      = command,
                             dependencies = dependencies)
