#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import glob
import random

import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions

from paleomix.node import CommandNode
from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder, \
    AtomicMPICmdBuilder, \
    use_customizable_cli_parameters, \
    create_customizable_cli_parameters

from paleomix.nodegraph import \
    FileStatusCache


EXAML_VERSION = versions.Requirement(call   = ("examl", "-version"),
                                     search = r"version (\d+)\.(\d+)\.(\d+)",
                                     checks = versions.GE(3, 0, 0))

PARSER_VERSION = versions.Requirement(call   = ("parse-examl", "-h"),
                                      search = r"version (\d+)\.(\d+)\.(\d+)",
                                      checks = versions.GE(3, 0, 0))


class ExaMLParserNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, input_partition, output_file, dependencies = ()):
        """
        Arguments:
        input_alignment  -- An alignment file in a format readable by RAxML.
        input_partition  -- A set of partitions in a format readable by RAxML.
        output_filename  -- Filename for the output binary sequence."""

        command = AtomicCmdBuilder("parse-examl", set_cwd = True)

        command.set_option("-s", "%(TEMP_OUT_ALN)s")
        command.set_option("-q", "%(TEMP_OUT_PART)s")
        # Output file will be named output.binary, and placed in the CWD
        command.set_option("-n", "output")

        # Substitution model
        command.set_option("-m", "DNA", fixed = False)


        command.set_kwargs(# Auto-delete: Symlinks
                          TEMP_OUT_PART   = os.path.basename(input_partition),
                          TEMP_OUT_ALN    = os.path.basename(input_alignment),

                          # Input files, are not used directly (see below)
                          IN_ALIGNMENT    = input_alignment,
                          IN_PARTITION    = input_partition,

                          # Final output file, are not created directly
                          OUT_BINARY      = output_file,

                          CHECK_EXAML     = PARSER_VERSION)

        return {"command" : command}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._symlinks = [os.path.abspath(parameters.input_alignment),
                          os.path.abspath(parameters.input_partition)]
        self._output_file = os.path.basename(parameters.output_file)


        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<ExaMLParser: '%s' -> '%s'>" \
                                 % (parameters.input_alignment, parameters.output_file),
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        # Required to avoid the creation of files outside the temp folder
        for filename in self._symlinks:
            source      = os.path.abspath(filename)
            destination = os.path.join(temp, os.path.basename(filename))

            os.symlink(source, destination)


    def _teardown(self, config, temp):
        os.remove(os.path.join(temp, "RAxML_info.output"))

        source      = os.path.join(temp, "output.binary")
        destination = fileutils.reroot_path(temp, self._output_file)
        fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)


class ExaMLNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_binary, initial_tree, output_template, threads = 1, dependencies = ()):
        """
        Arguments:
        input_binary  -- A binary alignment file in a format readable by ExaML.
        output_template  -- A template string used to construct final filenames. Should consist
                            of a full path, including a single '%s', which is replaced with the
                            variable part of RAxML output files (e.g. 'info', 'bestTree', ...).
                            Example destination: '/disk/project/SN013420.RAxML.%s'
                            Example output:      '/disk/project/SN013420.RAxML.bestTree'"""

        # TODO: Make MPIParams!
        command = AtomicMPICmdBuilder("examl", threads = threads)

        # Ensures that output is saved to the temporary directory
        command.set_option("-w", "%(TEMP_DIR)s")

        command.set_option("-s", "%(IN_ALN)s")
        command.set_option("-t", "%(IN_TREE)s")
        command.set_option("-n", "Pypeline")

        command.set_kwargs(IN_ALN=input_binary,
                           IN_TREE=initial_tree,

                           # Final output files, are not created directly
                           OUT_INFO=output_template % "info",
                           OUT_BESTTREE=output_template % "result",
                           OUT_BOOTSTRAP=output_template % "log",

                           # Only generated by newer versions of ExaML
                           TEMP_OUT_MODELFILE=os.path.basename(output_template
                                                               % "modelFile"),

                           CHECK_EXAML=EXAML_VERSION)

        # Use the GAMMA model of NT substitution by default
        command.set_option("-m", "GAMMA", fixed = False)

        return {"command"         : command}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._dirname  = os.path.dirname(parameters.output_template)
        self._template = os.path.basename(parameters.output_template)

        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<ExaML (%i thread(s)): '%s' -> '%s'>" \
                                 % (parameters.threads,
                                    parameters.input_binary,
                                    parameters.output_template),
                             threads      = parameters.threads,
                             dependencies = parameters.dependencies)

    def _create_temp_dir(self, _config):
        """Called by 'run' in order to create a temporary folder.
        To allow restarting from checkpoints, we use a fixed folder
        determined by the output_template."""
        temp = os.path.join(self._dirname, self._template % ("temp",))
        fileutils.make_dirs(temp)
        return temp

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        # The temp folder may contain old files:
        # Remove old pipes to prevent failure at _teardown
        for pipe_fname in glob.glob(os.path.join(temp, "pipe*")):
            fileutils.try_remove(pipe_fname)
        # ExaML refuses to overwrite old info files
        fileutils.try_remove(os.path.join(temp, "ExaML_info.Pypeline"))

        # Resume from last checkpoint, if one such was generated
        checkpoints = glob.glob(os.path.join(temp,
                                "ExaML_binaryCheckpoint.Pypeline_*"))
        if not checkpoints:
            return

        cache = FileStatusCache()
        if not cache.are_files_outdated(self.input_files, checkpoints):
            checkpoints.sort(key=lambda fname: int(fname.rsplit("_", 1)[-1]))

            # FIXME: Less hacky solution to modifying AtomicCmds needed
            self._command._command.append("-R")
            self._command._command.append(checkpoints[-1])
        else:
            for fpath in checkpoints:
                fileutils.try_remove(fpath)

    def _teardown(self, config, temp):
        for filename in os.listdir(temp):
            match = re.match("ExaML_(.*).Pypeline", filename)
            if match:
                if "binaryCheckpoint" in match.groups():
                    os.remove(os.path.join(temp, filename))
                else:
                    source      = os.path.join(temp, filename)
                    destination = os.path.join(temp, self._template % match.groups())

                    fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)




class ParsimonatorNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_alignment, output_tree, dependencies = ()):
        """
        Arguments:
        input_alignment  -- An alignment file in a format readable by RAxML.
        output_tree      -- Filename for the output newick tree."""

        command = AtomicCmdBuilder("parsimonator", set_cwd = True)

        command.set_option("-s", "%(TEMP_OUT_ALN)s")
        command.set_option("-n", "output")
        # Random seed for the stepwise addition process
        command.set_option("-p", int(random.random() * 2**31 - 1), fixed = False)

        command.set_kwargs(# Auto-delete: Symlinks
                          TEMP_OUT_ALN   = os.path.basename(input_alignment),

                          # Input files, are not used directly (see below)
                          IN_ALIGNMENT    = input_alignment,

                          # Final output file, are not created directly
                          OUT_TREE       = output_tree)

        return {"command"         : command}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._symlinks = [os.path.abspath(parameters.input_alignment)]
        self._output_tree = os.path.basename(parameters.output_tree)


        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<Parsimonator: '%s' -> '%s'>" \
                                 % (parameters.input_alignment, parameters.output_tree),
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        # Required to avoid the creation of files outside the temp folder
        for filename in self._symlinks:
            source      = os.path.abspath(filename)
            destination = os.path.join(temp, os.path.basename(filename))

            os.symlink(source, destination)


    def _teardown(self, config, temp):
        os.remove(os.path.join(temp, "RAxML_info.output"))

        source      = os.path.join(temp, "RAxML_parsimonyTree.output.0")
        destination = fileutils.reroot_path(temp, self._output_tree)
        fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)

