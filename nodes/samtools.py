import os

import node
from atomiccmd import *


class Genotype(node.Node):
    def __init__(self, config, destination, reference, infile, outfile, dependencies = []):
        cmd_pileup   = AtomicCmd(destination,
                                 ["samtools", "mpileup", "-uf", "%(IN_REFERENCE)s", "-EA", "%(IN_BAMFILE)s"],
                                 IN_REFERENCE = reference,
                                 IN_BAMFILE   = infile,
                                 stdout       = AtomicCmd.PIPE,
                                 stderr       = outfile + ".pileup_log")
        cmd_genotype = AtomicCmd(destination,
                                 ["bcftools", "view", "-bg", "-"],
                                 stdin        = cmd_pileup,
                                 stdout       = outfile,
                                 stderr       = outfile + ".genotype_log")

        description = "<Genotyper: '%s' -> '%s'>" \
            % (infile, os.path.join(destination, outfile))
            
        node.Node.__init__(self, 
                           description = description,
                           command = AtomicSet([cmd_pileup, cmd_genotype]),
                           dependencies = dependencies)
