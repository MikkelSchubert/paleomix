import os

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd



class JModelTestNode(CommandNode):
    pileup_args = "-EA"

    def __init__(self, config, infile, outfile, threads = 1, dependencies = ()):
        call = ["java", "-Djava.awt.headless=true",
                "-jar", os.path.join(config.jar_root, "jModelTest.jar"),
                "-f",         # Include models with unequals base frecuencies 
                "-i",         # Include models with a proportion invariable sites
                "-g", 8,      # Include models with rate variation among sites and number of categories
                "-s", 11,     # Number of substitution schemes
                "-S", "BEST", # Tree topology search operation option
                "-AICc",      # Calculate the corrected Akaike Information Criterion
                "-BIC",       # Calculate the Bayesian Information Criterion
                "-tr", threads,
		"-d", "%(IN_ALIGNMENT)s"]

        jmodeltest  = AtomicCmd(call,
                                IN_ALIGNMENT = infile,
                                OUT_STDOUT   = outfile)

        description = "<jModelTest: '%s' -> '%s'>" % (infile, outfile)
        CommandNode.__init__(self, 
                             description  = description,
                             command      = jmodeltest,
                             dependencies = dependencies)
