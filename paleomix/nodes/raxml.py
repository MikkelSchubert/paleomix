#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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

import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions

from paleomix.node import CommandNode
from paleomix.atomiccmd.builder import AtomicCmdBuilder


RAXML_VERSION = versions.Requirement(
    call=("raxmlHPC", "-version"),
    search=r"version (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(8, 2, 9),
)
RAXML_PTHREADS_VERSION = versions.Requirement(
    call=("raxmlHPC-PTHREADS", "-version"),
    search=r"version (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(8, 2, 9),
)


class RAxMLRapidBSNode(CommandNode):
    def __init__(
        self,
        input_alignment,
        output_template,
        input_partition=None,
        model="GTRGAMMAI",
        replicates="autoMRE",
        threads=1,
        dependencies=(),
    ):
        """
        Arguments:
        input_alignment  -- An alignment file in a format readable by RAxML.
        input_partition  -- A set of partitions in a format readable by RAxML.
        output_template  -- A template string used to construct final filenames. Should
                            consist of a full path, including a single '%s', which is
                            replaced with the  variable part of RAxML output files (e.g.
                            'info', 'bestTree', ...).

                            Example destination: '/disk/project/SN013420.RAxML.%s'
                            Example output:      '/disk/project/SN013420.RAxML.bestTree'
        """

        if threads > 1:
            command = AtomicCmdBuilder("raxmlHPC-PTHREADS")
            command.set_option("-T", threads)
            version = RAXML_PTHREADS_VERSION
        else:
            command = AtomicCmdBuilder("raxmlHPC")
            version = RAXML_VERSION

        # Perform rapid bootstrapping
        command.set_option("-f", "a")
        # Output files are saved with a .PALEOMIX postfix, and subsequently renamed
        command.set_option("-n", "PALEOMIX")
        # Ensures that output is saved to the temporary directory
        command.set_option("-w", "%(TEMP_DIR)s")
        # Symlink to sequence and partitions, to prevent the creation of *.reduced files
        # outside temp folder. In addition, it may be nessesary to remove the .reduced
        # files if created
        command.set_option("-s", "%(TEMP_OUT_ALN)s")

        if input_partition is not None:
            command.set_option("-q", "%(TEMP_OUT_PART)s")
            command.set_kwargs(
                IN_PARTITION=input_partition,
                TEMP_OUT_PART=os.path.basename(input_partition),
                TEMP_OUT_PART_R=os.path.basename(input_partition) + ".reduced",
            )

        command.set_kwargs(
            # Auto-delete: Symlinks and .reduced files that RAxML may generate
            TEMP_OUT_ALN=os.path.basename(input_alignment),
            TEMP_OUT_ALN_R=os.path.basename(input_alignment) + ".reduced",
            # Input files, are not used directly (see below)
            IN_ALIGNMENT=input_alignment,
            # Final output files, are not created directly
            OUT_INFO=output_template % "info",
            OUT_BESTTREE=output_template % "bestTree",
            OUT_BOOTSTRAP=output_template % "bootstrap",
            OUT_BIPART=output_template % "bipartitions",
            OUT_BIPARTLABEL=output_template % "bipartitionsBranchLabels",
            CHECK_VERSION=version,
        )

        # Use the GTRGAMMA model of NT substitution by default
        command.set_option("-m", model, fixed=False)
        # Enable Rapid Boostrapping and set random seed. May be set to a fixed value to
        # allow replicability.
        command.set_option("-x", int(random.random() * 2 ** 31 - 1), fixed=False)
        # Set random seed for parsimony inference. May be set to allow replicability.
        command.set_option("-p", int(random.random() * 2 ** 31 - 1), fixed=False)
        # Terminate bootstrapping upon convergence, not after N repetitions
        command.set_option("-N", replicates, fixed=False)

        self._symlinks = [input_alignment, input_partition]
        self._template = os.path.basename(output_template)

        CommandNode.__init__(
            self,
            command=command.finalize(),
            description="inferring phylogeny from %s using RAxML" % (input_alignment,),
            threads=threads,
            dependencies=dependencies,
        )

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        # Required to avoid the creation of files outside the temp folder
        for filename in self._symlinks:
            if filename is not None:
                source = os.path.abspath(filename)
                destination = os.path.join(temp, os.path.basename(filename))

                os.symlink(source, destination)

    def _teardown(self, config, temp):
        for filename in os.listdir(temp):
            match = re.match("RAxML_(.*).PALEOMIX", filename)
            if match:
                source = os.path.join(temp, filename)
                destination = os.path.join(temp, self._template % match.groups())

                fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)


class RAxMLParsimonyTreeNode(CommandNode):
    def __init__(self, input_alignment, input_partitions, output_tree, dependencies=()):
        command = AtomicCmdBuilder("raxmlHPC")

        # Compute a randomized parsimony starting tree
        command.set_option("-y")
        # Output files are saved with a .Pypeline postfix, and subsequently renamed
        command.set_option("-n", "Pypeline")
        # Model required, but not used
        command.set_option("-m", "GTRGAMMA")
        # Ensures that output is saved to the temporary directory
        command.set_option("-w", "%(TEMP_DIR)s")
        # Set random seed for bootstrap generation. May be set to allow replicability
        command.set_option("-p", int(random.random() * 2 ** 31 - 1), fixed=False)

        # Symlink to sequence and partitions, to prevent the creation of *.reduced files
        # outside temp folder
        command.set_option("-s", "%(TEMP_OUT_ALIGNMENT)s")
        command.set_option("-q", "%(TEMP_OUT_PARTITION)s")

        command.set_kwargs(
            IN_ALIGNMENT=input_alignment,
            IN_PARTITION=input_partitions,
            # TEMP_OUT_ is used to automatically remove these files
            TEMP_OUT_ALIGNMENT="RAxML_alignment",
            TEMP_OUT_PARTITION="RAxML_partitions",
            TEMP_OUT_INFO="RAxML_info.Pypeline",
            OUT_TREE=output_tree,
            CHECK_VERSION=RAXML_VERSION,
        )

        self._input_alignment = input_alignment
        self._input_partitions = input_partitions
        self._output_tree = output_tree

        CommandNode.__init__(
            self,
            command=command.finalize(),
            description="inferring parsimony phylogeny from %s using RAxML"
            % (input_alignment,),
            dependencies=dependencies,
        )

    def _setup(self, config, temp):
        os.symlink(
            os.path.abspath(self._input_alignment),
            os.path.join(temp, "RAxML_alignment"),
        )
        os.symlink(
            os.path.abspath(self._input_partitions),
            os.path.join(temp, "RAxML_partitions"),
        )
        CommandNode._setup(self, config, temp)

    def _teardown(self, config, temp):
        basename = os.path.basename(self._output_tree)
        os.rename(
            os.path.join(temp, "RAxML_parsimonyTree.Pypeline"),
            os.path.join(temp, basename),
        )

        CommandNode._teardown(self, config, temp)
