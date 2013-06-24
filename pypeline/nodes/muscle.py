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

from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd.command import AtomicCmd


class MuscleNode(CommandNode):
    def __init__(self, infile, prefix, dependencies = ()):
        call = ["muscle",
                "-quiet",
                "-maxiters", 32,
                "-in", "%(IN_FASTA)s",
                "-fastaout", "%(OUT_AFA)s",
                "-phyiout",  "%(OUT_PHYS)s"]

        command = AtomicCmd(call,
                            IN_FASTA = infile,
                            OUT_AFA  = prefix + ".afa",
                            OUT_PHYS = prefix + ".phy")

        description = "<MuscleNode: '%s' -> '%s'>" \
                % (infile, prefix + ".*")

        CommandNode.__init__(self,
                             command      = command,
                             description  = description,
                             dependencies = dependencies)
 



class MetaMuscleNode(MetaNode):
    def __init__(self, rootdir, sequences, dependencies = ()):
        subnodes = []
        for sequence in sequences:
            prefix  = os.path.join(rootdir, sequence)
            node    = MuscleNode(infile       = prefix + ".fasta",
                                 prefix       = prefix, 
                                 dependencies = dependencies)

            subnodes.append(node)

        MetaNode.__init__(self,
                          description  = "<MuscleAlignSequences: In '%s'>" \
                              % (rootdir,),
                          subnodes     = subnodes,
                          dependencies = dependencies)
