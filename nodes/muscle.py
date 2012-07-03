import os

from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd import AtomicCmd


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
                          description = "<MuscleAlignSequences: In '%s'>" \
                                % (rootdir,),
                          subnodes    = subnodes)                            
