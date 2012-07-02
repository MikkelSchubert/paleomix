import os

import pypeline.fileutils as fileutils

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd


class MuscleNode(CommandNode):
    def __init__(self, config, infile, outfile, dependencies = ()):
        call = ["muscle", 
                "-maxiters", 32,
                "-in", "%(IN_FASTA)s",
                "-physout", "%(OUT_PHYS)s"]

        command = AtomicCmd(call,
                            IN_FASTA = infile,
                            OUT_PHYS = outfile)

        description = "<MuscleNode: '%s' -> '%s'>" \
                % (infile, outfile)

        CommandNode.__init__(self,
                             command      = command,
                             description  = description,
                             dependencies = dependencies)
                             
