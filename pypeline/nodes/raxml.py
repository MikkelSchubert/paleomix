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
import re
import random

import pypeline.common.fileutils as fileutils

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd




class RAxMLReduceNode(CommandNode):
    def __init__(self, inalignment, outalignment, inpartitions, outpartitions, dependencies = ()):
        call = ["raxmlHPC",
                "-f", "c",
                "-m", "GTRGAMMA", # Model required, but not used
                "-n", "GTRGAMMA",
                "-w", "%(TEMP_DIR)s",
                "-s", "%(TEMP_IN_ALIGNMENT)s",
                "-q", "%(TEMP_IN_PARTITIONS)s"]

        self._kwargs = {"IN_ALIGNMENT"        : inalignment,
                        "IN_PARTITIONS"       : inpartitions,
                            
                        "TEMP_IN_ALIGNMENT"   : "RAXML_alignment",
                        "TEMP_IN_PARTITIONS"  : "RAXML_partitions",
                        "TEMP_OUT_INFO"       : "RAxML_info.GTRGAMMA",

                        "OUT_ALIGNMENT"       : outalignment,
                        "OUT_PARTITIONS"      : outpartitions}

        CommandNode.__init__(self,
                             command      = AtomicCmd(call, **self._kwargs),
                             description  = "<RAxMLReduce: '%s' -> '%s'>" \
                                     % (inalignment, outalignment),
                             dependencies = dependencies)


    def _setup(self, config, temp):
        for key in ("IN_ALIGNMENT", "IN_PARTITIONS"):
            source      = self._kwargs[key]
            destination = os.path.join(temp, self._kwargs["TEMP_" + key])

            fileutils.copy_file(source, destination)

        CommandNode._setup(self, config, temp)


    def _teardown(self, config, temp):
        for postfix in ("ALIGNMENT", "PARTITIONS"):
            filenames = [self._kwargs["TEMP_IN_" + postfix],
                         self._kwargs["TEMP_IN_" + postfix] + ".reduced",
                         self._kwargs["OUT_" + postfix]]

            for (source, destination) in zip(filenames, filenames[1:]):
                source      = fileutils.reroot_path(temp, source)
                destination = fileutils.reroot_path(temp, destination)

                if not os.path.exists(destination):
                    fileutils.move_file(source, destination)
                elif source != destination:
                    os.remove(source)
        
        CommandNode._teardown(self, config, temp)




class RAxMLRapidBSNode(CommandNode):
    def __init__(self, infile, partitions, destination, model = "GTRGAMMAI", replicates = "autoMRE", bs_branchln = False, threads = 1, dependencies = ()):
        """ 
        Arguments:
        infile      -- An alignment file in a format readable by RAxML.
        partitions  -- A set of partitions in a format readable by RAxML.
        destination -- A template string used to construct final filenames. Should consist
                       of a full path, including a single '%s', which is replaced with the
                       variable part of RAxML output files (e.g. 'info', 'bestTree', ...).
                       Example destination: '/disk/project/SN013420.RAxML.%s'
                       Example output:      '/disk/project/SN013420.RAxML.bestTree'
        model       -- DNA or Amino acid substitution model to use.
        replicates  -- Number of bootstrap replicates. Defaults to automatic selection.
        bs_branchln -- Calculate branch lengths for bootstrap trees.
        threads     -- Number of threads to use for each RAxML instance."""

        self._symlinks = [infile, partitions]
        self._template = os.path.basename(destination)

        call = ["raxmlHPC-PTHREADS" if (threads > 1) else "raxmlHPC",
                "-f", "a",
                "-m", model,
                "-n", "RapidBS",
                "-s", "%(TEMP_OUT_ALN)s",
                "-q", "%(TEMP_OUT_PART)s",
                "-w", "%(TEMP_DIR)s",
                "-x", int(random.random() * 2**31 - 1),
                "-p", int(random.random() * 2**31 - 1),
                "-N", replicates]
        if bs_branchln:
            call.append("-k")
        if threads > 1:
            call.extend(("-T", threads))


        command = AtomicCmd(call,
                            IN_ALIGNMENT    = infile,
                            IN_PARTITIONS   = partitions,

                            OUT_INFO        = destination % "info",
                            OUT_BESTTREE    = destination % "bestTree",
                            OUT_BOOTSTRAP   = destination % "bootstrap",
                            OUT_BIPART      = destination % "bipartitions",
                            OUT_BIPARTLABEL = destination % "bipartitionsBranchLabels",

                            # RAxML may write reduced alignment. These are not saved. If
                            # needed, they may be generated upfront using RAxMLReduceNode.
                            TEMP_OUT_ALN    = os.path.basename(infile),
                            TEMP_OUT_PART   = os.path.basename(partitions),
                            TEMP_OUT_R_ALN  = os.path.basename(infile) + ".reduced",
                            TEMP_OUT_R_PART = os.path.basename(partitions) + ".reduced")

        CommandNode.__init__(self,
                             command      = command,
                             description  = "<RAxMLRapidBS (%i threads): '%s' -> '%s'>" \
                                     % (threads, infile, destination),
                             threads      = threads,
                             dependencies = dependencies)

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        # Required to avoid the creation of files outside the temp folder
        for filename in self._symlinks:
            source      = os.path.abspath(filename)
            destination = os.path.join(temp, os.path.basename(filename)) 

            os.symlink(source, destination)

    
    def _teardown(self, config, temp):
        for filename in os.listdir(temp):
            match = re.match("RAxML_(.*).RapidBS", filename)
            if match:
                source      = os.path.join(temp, filename)
                destination = os.path.join(temp, self._template % match.groups())

                fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)
