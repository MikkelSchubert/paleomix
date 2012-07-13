import os

from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd import AtomicCmd


# Presets mainly taken from
# http://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html
_PRESETS = {
    "auto"     : ["mafft", "--auto"],
    "FFT-NS-1" : ["fftns", "--retree", 1],
    "FFT-NS-2" : ["fftns"],
    "FFT-NS-i" : ["fftnsi", "--maxiterate", 1000],
    "NW-NS-i"  : ["nwnsi",  "--maxiterate", 1000],
    "L-INS-i"  : ["linsi",  "--maxiterate", 1000],
    "E-INS-i"  : ["einsi",  "--maxiterate", 1000],
    "G-INS-i"  : ["ginsi",  "--maxiterate", 1000]}
    


class MAFFTNode(CommandNode):
    def __init__(self, infile, outfile, preset = "auto", dependencies = ()):
        call = _PRESETS[preset] + ["--quiet", "%(IN_FASTA)s"]
        command = AtomicCmd(call,
                            IN_FASTA   = infile,
                            OUT_STDOUT = outfile)

        description = "<MAFFTNode: '%s' -> '%s'>" \
                % (infile, outfile)

        CommandNode.__init__(self,
                             command      = command,
                             description  = description,
                             dependencies = dependencies)
 



class MetaMAFFTNode(MetaNode):
    def __init__(self, rootdir, sequences, preset = "auto", dependencies = ()):
        subnodes = []
        for sequence in sequences:
            prefix  = os.path.join(rootdir, sequence)
            node    = MAFFTNode(infile       = prefix + ".fasta",
                                outfile      = prefix + ".afa",
                                preset       = preset,
                                dependencies = dependencies)

            subnodes.append(node)

        MetaNode.__init__(self,
                          description = "<MAFFTAlignSequences: In '%s'>" \
                                % (rootdir,),
                          subnodes    = subnodes)                            
