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
import os

import pypeline.common.fileutils as fileutils
from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.common.fileutils import swap_ext


class _IndelTrainerNode(CommandNode):
    def __init__(self, config, reference, infile, outfile, dependencies = ()):
        call  = ["java", "-jar", config.gatk_jar,
                 "-T", "RealignerTargetCreator", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-o", "%(OUT_INTERVALS)s"]

        command = AtomicCmd(call, 
                            IN_REFERENCE  = reference,
                            IN_BAMFILE    = infile,
                            IN_BAIFILE    = swap_ext(infile, ".bai"),
                            OUT_INTERVALS = outfile)
        
        description = "<Train Indel Realigner: '%s' -> '%s'>" \
            % (infile, outfile)

        CommandNode.__init__(self, 
                             description = description,
                             command = command,
                             dependencies = dependencies)



class _IndelRealignerNode(CommandNode):
    def __init__(self, config, reference, intervals, infile, outfile, dependencies = ()):
        call  = ["java", "-jar", config.gatk_jar,
                 "-T", "IndelRealigner", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-targetIntervals", "%(IN_INTERVALS)s",
                 "-o", "%(OUT_BAMFILE)s"]
        
        command = AtomicCmd(call, 
                            IN_REFERENCE = reference,
                            IN_BAMFILE   = infile,
                            IN_BAIFILE   = swap_ext(infile, ".bai"),
                            IN_INTERVALS = intervals,
                            OUT_BAMFILE  = outfile,
                            OUT_INDEX    = fileutils.swap_ext(outfile, ".bai"))

        description = "<Indel Realign: '%s' -> '%s'>" \
            % (infile, outfile)

        CommandNode.__init__(self, 
                             description = description,
                             command = command,
                             dependencies = dependencies)



class IndelRealignerNode(MetaNode):
    def __init__(self, config, reference, infile, outfile, intervals = None, dependencies = ()):
        if not intervals:
            intervals = outfile + ".intervals"

        trainer = _IndelTrainerNode(config         = config,
                                    reference      = reference, 
                                    infile         = infile,
                                    outfile        = intervals,
                                    dependencies   = dependencies)
        aligner = _IndelRealignerNode(config       = config,
                                      reference    = reference, 
                                      intervals    = intervals,
                                      infile       = infile, 
                                      outfile      = outfile,
                                      dependencies = trainer)
        
        MetaNode.__init__(self, 
                          description  = "<GATK Indel Realigner: '%s' -> '%s'>" \
                              % (infile, outfile),
                          subnodes     = [trainer, aligner],
                          dependencies = dependencies)
