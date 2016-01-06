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

import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions

from paleomix.node import CommandNode
from paleomix.atomiccmd.builder import \
     AtomicCmdBuilder, \
     use_customizable_cli_parameters, \
     create_customizable_cli_parameters


RAXML_VERSION = versions.Requirement(call   = ("raxmlHPC", "-version"),
                                     search = r"version (\d+)\.(\d+)\.(\d+)",
                                     checks = versions.GE(7, 3, 2))
RAXML_PTHREADS_VERSION = versions.Requirement(call   = ("raxmlHPC-PTHREADS", "-version"),
                                              search = r"version (\d+)\.(\d+)\.(\d+)",
                                              checks = versions.GE(7, 3, 2))


class RAxMLReduceNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, input_partition, output_alignment, output_partition, dependencies = ()):
        command = AtomicCmdBuilder("raxmlHPC")

        # Read and (in the case of empty columns) reduce input
        command.set_option("-f", "c")
        # Output files are saved with a .Pypeline postfix, and subsequently renamed
        command.set_option("-n", "Pypeline")
        # Model required, but not used
        command.set_option("-m", "GTRGAMMA")
        # Ensures that output is saved to the temporary directory
        command.set_option("-w", "%(TEMP_DIR)s")

        # Symlink to sequence and partitions, to prevent the creation of *.reduced files outside temp folder
        # In addition, it may be nessesary to remove the .reduced files if created
        command.set_option("-s", "%(TEMP_IN_ALIGNMENT)s")
        command.set_option("-q", "%(TEMP_IN_PARTITION)s")

        command.set_kwargs(IN_ALIGNMENT      = input_alignment,
                          IN_PARTITION      = input_partition,

                          TEMP_IN_ALIGNMENT = "RAxML_alignment",
                          TEMP_IN_PARTITION = "RAxML_partitions",
                          TEMP_OUT_INFO     = "RAxML_info.Pypeline",

                          OUT_ALIGNMENT     = output_alignment,
                          OUT_PARTITION     = output_partition,
                          CHECK_VERSION     = RAXML_VERSION)

        return {"command" : command}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._kwargs = parameters.command.kwargs
        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<RAxMLReduce: '%s' -> '%s'>" \
                                     % (parameters.input_alignment, parameters.output_alignment),
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        for key in ("IN_ALIGNMENT", "IN_PARTITION"):
            source      = os.path.abspath(self._kwargs[key])
            destination = os.path.join(temp, self._kwargs["TEMP_" + key])

            os.symlink(source, destination)

        CommandNode._setup(self, config, temp)


    def _teardown(self, config, temp):
        for postfix in ("ALIGNMENT", "PARTITION"):
            filenames = [self._kwargs["TEMP_IN_" + postfix],
                         self._kwargs["TEMP_IN_" + postfix] + ".reduced",
                         self._kwargs["OUT_" + postfix]]

            for (source, destination) in zip(filenames, filenames[1:]):
                source      = fileutils.reroot_path(temp, source)
                destination = fileutils.reroot_path(temp, destination)

                if not os.path.exists(destination):
                    fileutils.copy_file(source, destination)
                os.remove(source)

        CommandNode._teardown(self, config, temp)


class RAxMLBootstrapNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, input_partition, template, start = 0, bootstraps = 50, dependencies = ()):
        command = AtomicCmdBuilder("raxmlHPC", set_cwd = True)

        # Read and (in the case of empty columns) reduce input
        command.set_option("-f", "j")
        # Output files are saved with a .Pypeline postfix, and subsequently renamed
        command.set_option("-n", "Pypeline")
        # Model required, but not used
        command.set_option("-m", "GTRGAMMA")
        # Set random seed for bootstrap generation. May be set to a fixed value to allow replicability.
        command.set_option("-b", int(random.random() * 2**31 - 1), fixed = False)
        # Generate a single bootstrap alignment (makes growing the number of bootstraps easier).
        command.set_option("-N", int(bootstraps), fixed = False)

        # Symlink to sequence and partitions, to prevent the creation of *.reduced files outside temp folder
        # In addition, it may be nessesary to remove the .reduced files if created
        command.set_option("-s", "input.alignment")
        command.set_option("-q", "input.partition")

        bootstrap_files = {"IN_ALIGNMENT" : input_alignment,
                           "IN_PARTITION" : input_partition,
                           "TEMP_OUT_INF" : "RAxML_info.Pypeline",
                           "TEMP_OUT_ALN" : "input.alignment",
                           "TEMP_OUT_PAR" : "input.partition",
                           "CHECK_VERSION": RAXML_VERSION}

        for (index, (_, filename)) in enumerate(cls._bootstraps(template, bootstraps, start)):
            bootstrap_files["OUT_BS_%03i" % index] = filename
        command.set_kwargs(**bootstrap_files)

        return {"command" : command}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._input_alignment = parameters.input_alignment
        self._input_partition = parameters.input_partition
        self._output_template = parameters.template
        self._bootstrap_num   = parameters.bootstraps
        self._bootstrap_start = parameters.start

        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<RAxMLBootstrap: '%s' -> '%s' (%i .. %i>" \
                                     % (parameters.input_alignment, parameters.template,
                                        parameters.start, parameters.start + parameters.bootstraps - 1),
                             dependencies = parameters.dependencies)

    def _setup(self, config, temp):
        os.symlink(os.path.realpath(self._input_alignment), os.path.join(temp, "input.alignment"))
        os.symlink(os.path.realpath(self._input_partition), os.path.join(temp, "input.partition"))


    def _teardown(self, config, temp):
        template   = self._output_template
        bootstraps = self._bootstrap_num
        start      = self._bootstrap_start
        for (src_file, dst_file) in self._bootstraps(template, bootstraps, start):
            src_file = os.path.join(temp, src_file)
            dst_file = fileutils.reroot_path(temp, dst_file)
            fileutils.move_file(src_file, dst_file)
        CommandNode._teardown(self, config, temp)

    @classmethod
    def _bootstraps(cls, template, number, start):
        for bootstrap in range(number):
            src_file = "input.alignment.BS%i" % (bootstrap,)
            dst_file = template % (bootstrap + start,)
            yield (src_file, dst_file)


class RAxMLRapidBSNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, output_template, input_partition=None,
                  threads=1, dependencies=()):
        """
        Arguments:
        input_alignment  -- An alignment file in a format readable by RAxML.
        input_partition  -- A set of partitions in a format readable by RAxML.
        output_template  -- A template string used to construct final filenames. Should consist
                            of a full path, including a single '%s', which is replaced with the
                            variable part of RAxML output files (e.g. 'info', 'bestTree', ...).
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
        # Symlink to sequence and partitions, to prevent the creation of *.reduced files outside temp folder
        # In addition, it may be nessesary to remove the .reduced files if created
        command.set_option("-s", "%(TEMP_OUT_ALN)s")

        if input_partition is not None:
            command.set_option("-q", "%(TEMP_OUT_PART)s")
            command.set_kwargs(IN_PARTITION=input_partition,
                               TEMP_OUT_PART=os.path.basename(input_partition),
                               TEMP_OUT_PART_R=os.path.basename(input_partition) + ".reduced")

        command.set_kwargs(  # Auto-delete: Symlinks and .reduced files that RAxML may generate
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

                           CHECK_VERSION=version)

        # Use the GTRGAMMA model of NT substitution by default
        command.set_option("-m", "GTRGAMMAI", fixed=False)
        # Enable Rapid Boostrapping and set random seed. May be set to a fixed value to allow replicability.
        command.set_option("-x", int(random.random() * 2**31 - 1), fixed=False)
        # Set random seed for parsimony inference. May be set to a fixed value to allow replicability.
        command.set_option("-p", int(random.random() * 2**31 - 1), fixed=False)
        # Terminate bootstrapping upon convergence, rather than after a fixed number of repetitions
        command.set_option("-N", "autoMRE", fixed=False)

        return {"command": command}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._symlinks = [parameters.input_alignment,
                          parameters.input_partition]
        self._template = os.path.basename(parameters.output_template)

        CommandNode.__init__(self,
                             command=parameters.command.finalize(),
                             description="<RAxMLRapidBS: '%s' -> '%s'>"
                             % (parameters.input_alignment,
                                parameters.output_template % ("*",)),
                             threads=parameters.threads,
                             dependencies=parameters.dependencies)

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
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, input_partitions, output_tree, dependencies = ()):
        command = AtomicCmdBuilder("raxmlHPC")

        # Compute a randomized parsimony starting tree
        command.set_option("-y")
        # Output files are saved with a .Pypeline postfix, and subsequently renamed
        command.set_option("-n", "Pypeline")
        # Model required, but not used
        command.set_option("-m", "GTRGAMMA")
        # Ensures that output is saved to the temporary directory
        command.set_option("-w", "%(TEMP_DIR)s")
        # Set random seed for bootstrap generation. May be set to a fixed value to allow replicability.
        command.set_option("-p", int(random.random() * 2**31 - 1), fixed = False)

        # Symlink to sequence and partitions, to prevent the creation of *.reduced files outside temp folder
        command.set_option("-s", "%(TEMP_OUT_ALIGNMENT)s")
        command.set_option("-q", "%(TEMP_OUT_PARTITION)s")

        command.set_kwargs(IN_ALIGNMENT       = input_alignment,
                           IN_PARTITION       = input_partitions,

                           # TEMP_OUT_ is used to automatically remove these files
                           TEMP_OUT_ALIGNMENT = "RAxML_alignment",
                           TEMP_OUT_PARTITION = "RAxML_partitions",
                           TEMP_OUT_INFO      = "RAxML_info.Pypeline",

                           OUT_TREE           = output_tree,

                           CHECK_VERSION      = RAXML_VERSION)

        return {"command" : command}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._input_alignment  = parameters.input_alignment
        self._input_partitions = parameters.input_partitions
        self._output_tree      = parameters.output_tree

        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<RAxMLParsimonyTree: '%s' -> '%s'>" \
                                     % (parameters.input_alignment, parameters.output_tree),
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        os.symlink(os.path.abspath(self._input_alignment),  os.path.join(temp, "RAxML_alignment"))
        os.symlink(os.path.abspath(self._input_partitions), os.path.join(temp, "RAxML_partitions"))
        CommandNode._setup(self, config, temp)


    def _teardown(self, config, temp):
        basename = os.path.basename(self._output_tree)
        os.rename(os.path.join(temp, "RAxML_parsimonyTree.Pypeline"),
                  os.path.join(temp, basename))

        CommandNode._teardown(self, config, temp)
