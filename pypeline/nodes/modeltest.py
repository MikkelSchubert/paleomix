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

from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd



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
