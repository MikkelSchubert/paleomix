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
from pypeline.atomicparams import *




class RAxMLReduceNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, input_partition, output_alignment, output_partition, dependencies = ()):
        command = AtomicParams("raxmlHPC")

        # Read and (in the case of empty columns) reduce input
        command.set_parameter("-f", "c")
        # Output files are saved with a .Pypeline postfix, and subsequently renamed
        command.set_parameter("-n", "Pypeline")
        # Model required, but not used
        command.set_parameter("-m", "GTRGAMMA")
        # Ensures that output is saved to the temporary directory
        command.set_parameter("-w", "%(TEMP_DIR)s")
        
        # Symlink to sequence and partitions, to prevent the creation of *.reduced files outside temp folder
        # In addition, it may be nessesary to remove the .reduced files if created
        command.set_parameter("-s", "%(TEMP_IN_ALIGNMENT)s")
        command.set_parameter("-q", "%(TEMP_IN_PARTITION)s")

        command.set_paths(IN_ALIGNMENT      = input_alignment,
                          IN_PARTITION      = input_partition,
                          
                          TEMP_IN_ALIGNMENT = "RAXML_alignment",
                          TEMP_IN_PARTITION = "RAXML_partitions",
                          TEMP_OUT_INFO     = "RAxML_info.GTRGAMMA",
                          
                          OUT_ALIGNMENT     = output_alignment,
                          OUT_PARTITION     = output_partition)

        return {"command" : command} 


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._kwargs = parameters.command.paths
        CommandNode.__init__(self,
                             command      = parameters.command.create_cmd(),
                             description  = "<RAxMLReduce: '%s' -> '%s'>" \
                                     % (parameters.input_alignment, parameters.output_alignment),
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        for key in ("IN_ALIGNMENT", "IN_PARTITION"):
            source      = self._kwargs[key]
            destination = os.path.join(temp, self._kwargs["TEMP_" + key])

            fileutils.copy_file(source, destination)

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
                    fileutils.move_file(source, destination)
                elif source != destination:
                    os.remove(source)
        
        CommandNode._teardown(self, config, temp)




class RAxMLRapidBSNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, input_partition, output_template, threads = 1, dependencies = ()):
        """ 
        Arguments:
        input_alignment  -- An alignment file in a format readable by RAxML.
        input_partition  -- A set of partitions in a format readable by RAxML.
        output_template  -- A template string used to construct final filenames. Should consist
                            of a full path, including a single '%s', which is replaced with the
                            variable part of RAxML output files (e.g. 'info', 'bestTree', ...).
                            Example destination: '/disk/project/SN013420.RAxML.%s'
                            Example output:      '/disk/project/SN013420.RAxML.bestTree'"""

        if threads > 1:
            command = AtomicParams("raxmlHPC-PTHREADS")
            command.set_parameter("-T", threads)
        else:
            command = AtomicParams("raxmlHPC")
        
        # Perform rapid bootstrapping
        command.set_parameter("-f", "a")
        # Output files are saved with a .Pypeline postfix, and subsequently renamed
        command.set_parameter("-n", "Pypeline")
        # Ensures that output is saved to the temporary directory
        command.set_parameter("-w", "%(TEMP_DIR)s")
        # Symlink to sequence and partitions, to prevent the creation of *.reduced files outside temp folder
        # In addition, it may be nessesary to remove the .reduced files if created
        command.set_parameter("-s", "%(TEMP_OUT_ALN)s")
        command.set_parameter("-q", "%(TEMP_OUT_PART)s")

        command.set_paths(# Auto-delete: Symlinks and .reduced files that RAxML may generate
                          TEMP_OUT_PART   = os.path.basename(input_partition),
                          TEMP_OUT_PART_R = os.path.basename(input_partition) + ".reduced",
                          TEMP_OUT_ALN    = os.path.basename(input_alignment),
                          TEMP_OUT_ALN_R  = os.path.basename(input_alignment) + ".reduced",
        
                          # Input files, are not used directly (see below)
                          IN_ALIGNMENT    = input_alignment,
                          IN_PARTITION    = input_partition,

                          # Final output files, are not created directly
                          OUT_INFO        = output_template % "info",
                          OUT_BESTTREE    = output_template % "bestTree",
                          OUT_BOOTSTRAP   = output_template % "bootstrap",
                          OUT_BIPART      = output_template % "bipartitions",
                          OUT_BIPARTLABEL = output_template % "bipartitionsBranchLabels")

        # Use the GTRGAMMAI model of NT substitution by default
        command.set_parameter("-m", "GTRGAMMAI", fixed = False)
        # Enable Rapid Boostrapping and set random seed. May be set to a fixed value to allow replicability.
        command.set_parameter("-x", int(random.random() * 2**31 - 1), fixed = False)
        # Set random seed for parsimony inference. May be set to a fixed value to allow replicability.
        command.set_parameter("-p", int(random.random() * 2**31 - 1), fixed = False)
        # Terminate bootstrapping upon convergence, rather than after a fixed number of repetitions      
        command.set_parameter("-N", "autoMRE", fixed = False) 
        
        return {"command"         : command}
    

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._symlinks = [parameters.input_alignment, 
                          parameters.input_partition]
        self._template = os.path.basename(parameters.output_template)


        CommandNode.__init__(self,
                             command      = parameters.command.create_cmd(),
                             description  = "<RAxMLRapidBS: '%s' -> '%s'>" \
                                 % (parameters.input_alignment, parameters.output_template),
                             threads      = parameters.threads,
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        # Required to avoid the creation of files outside the temp folder
        for filename in self._symlinks:
            source      = os.path.abspath(filename)
            destination = os.path.join(temp, os.path.basename(filename)) 

            os.symlink(source, destination)

    
    def _teardown(self, config, temp):
        for filename in os.listdir(temp):
            match = re.match("RAxML_(.*).Pypeline", filename)
            if match:
                source      = os.path.join(temp, filename)
                destination = os.path.join(temp, self._template % match.groups())

                fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)
