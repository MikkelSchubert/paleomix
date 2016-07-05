#!/usr/bin/python
#
# Copyright (c) 2016 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import os

import paleomix.common.rtools as rtools
import paleomix.tools.factory as factory

from paleomix.atomiccmd.command import AtomicCmd
from paleomix.common.formats.newick import Newick
from paleomix.node import CommandNode

from paleomix.tools.zonkey.common import RSCRIPT_VERSION


class MitoConsensusNode(CommandNode):
    def __init__(self, database, bamfile, output_prefix, dependencies=()):
        cmd = factory.new("zonkey_mito")
        cmd.add_value("%(IN_DATABASE)s")
        cmd.add_value("%(IN_BAMFILE)s")
        cmd.add_value("%(TEMP_OUT_PREFIX)s")

        cmd.set_kwargs(IN_DATABASE=database,
                       IN_BAMFILE=bamfile,
                       TEMP_OUT_PREFIX=os.path.basename(output_prefix),
                       OUT_PHYLIP=output_prefix + ".phy",
                       OUT_FASTA=output_prefix + ".fasta",
                       OUT_SUMMARY=output_prefix + ".summary")

        CommandNode.__init__(self,
                             description="<MitoConsensus -> '%s.*'>"
                             % (output_prefix,),
                             command=cmd.finalize(),
                             dependencies=dependencies)


class DrawPhylogenyNode(CommandNode):
    def __init__(self, samples, treefile, bootstraps, output_prefix,
                 dependencies=()):
        rscript = rtools.rscript("zonkey", "tinytree.r")

        cmd = AtomicCmd(("Rscript", rscript,
                         "%(TEMP_OUT_FILE)s",
                         "%(IN_SAMPLES)s",
                         "%(TEMP_OUT_PREFIX)s"),
                        AUX_RSCRIPT=rscript,
                        IN_SAMPLES=samples,
                        IN_FILE=treefile,
                        IN_BOOTSTRAPS=bootstraps,
                        TEMP_OUT_FILE="rerooted.newick",
                        TEMP_OUT_PREFIX=os.path.basename(output_prefix),
                        OUT_TREE_PDF=output_prefix + ".pdf",
                        OUT_TREE_PNG=output_prefix + ".png",
                        CHECK_RSCRIPT=RSCRIPT_VERSION,
                        CHECK_RSCRIPT_APE=rtools.requirement("ape"),
                        CHECK_RSCRIPT_GGPLOT2=rtools.requirement("ggplot2"),
                        CHECK_RSCRIPT_GRID=rtools.requirement("grid"))

        self._treefile = treefile
        self._bootstraps = bootstraps

        CommandNode.__init__(self,
                             description="<DrawPhylogeny -> '%s.*'>"
                             % (output_prefix,),
                             command=cmd,
                             dependencies=dependencies)

    def _setup(self, config, temp):
        with open(self._bootstraps) as handle:
            bootstraps = [Newick.from_string(line.strip())
                          for line in handle]

        with open(self._treefile) as handle:
            tree = Newick.from_string(handle.read().strip())

        tree = tree.reroot_on_midpoint()
        tree = tree.add_support(bootstraps, "{Percentage:.0f}")
        with open(os.path.join(temp, "rerooted.newick"), "w") as handle:
            handle.write("{}\n".format(tree))

        CommandNode._setup(self, config, temp)
