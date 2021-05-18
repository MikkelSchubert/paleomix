#!/usr/bin/python3
#
# Copyright (c) 2013 Mikkel Schubert <MikkelSch@gmail.com>
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

import paleomix.common.fileutils as fileutils
import paleomix.common.versions as versions

from paleomix.node import CommandNode
from paleomix.atomiccmd.command import AtomicCmd, InputFile, OutputFile, Executable

from paleomix.nodegraph import FileStatusCache


EXAML_VERSION = versions.Requirement(
    call=("examl", "-version"),
    search=r"version (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(3, 0, 21),
)

PARSER_VERSION = versions.Requirement(
    call=("parse-examl", "-h"),
    search=r"version (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(3, 0, 21),
)


class ExaMLParserNode(CommandNode):
    def __init__(self, input_alignment, input_partition, output_file, dependencies=()):
        """
        Arguments:
        input_alignment  -- An alignment file in a format readable by RAxML.
        input_partition  -- A set of partitions in a format readable by RAxML.
        output_filename  -- Filename for the output binary sequence."""

        command = AtomicCmd(
            "parse-examl",
            extra_files=[
                # Input files, are not used directly (see below)
                InputFile(input_alignment),
                InputFile(input_partition),
                # Final output file, are not created directly
                OutputFile(output_file),
            ],
            set_cwd=True,
            requirements=[PARSER_VERSION],
        )

        command.append_options(
            {
                # Auto-delete symlinks
                "-s": OutputFile(os.path.basename(input_alignment), temporary=True),
                "-q": OutputFile(os.path.basename(input_partition), temporary=True),
                # Output file will be named output.binary, and placed in the CWD
                "-n": "output",
                # Substitution model
                "-m": "DNA",
            }
        )

        CommandNode.__init__(
            self,
            command=command,
            description="pre-parsing %s for ExaML" % input_alignment,
            dependencies=dependencies,
        )

        self._symlinks = [
            os.path.abspath(input_alignment),
            os.path.abspath(input_partition),
        ]
        self._output_file = os.path.basename(output_file)

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        # Required to avoid the creation of files outside the temp folder
        for filename in self._symlinks:
            source = os.path.abspath(filename)
            destination = os.path.join(temp, os.path.basename(filename))

            os.symlink(source, destination)

    def _teardown(self, config, temp):
        os.remove(os.path.join(temp, "RAxML_info.output"))

        source = os.path.join(temp, "output.binary")
        destination = fileutils.reroot_path(temp, self._output_file)
        fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)


class ExaMLNode(CommandNode):
    def __init__(
        self,
        input_binary,
        initial_tree,
        output_template,
        model="GAMMA",
        threads=1,
        dependencies=(),
    ):
        """
        Arguments:
        input_binary  -- A binary alignment file in a format readable by ExaML.
        output_template  -- A template string used to construct final filenames. Should
                            consist of a full path, including a single '%s', which is
                            replaced with the variable part of RAxML output files (e.g.
                            'info', 'bestTree', ...).

                            Example destination: '/disk/project/SN013420.RAxML.%s'
                            Example output:      '/disk/project/SN013420.RAxML.bestTree'
        """
        if threads > 1:
            call = ["mpirun", "-n", threads, Executable("examl")]
        elif threads == 1:
            call = ["examl"]
        else:
            raise ValueError("threads must be 1 or greater, not %i" % threads)

        command = AtomicCmd(
            call,
            extra_files=[
                # Final output files, are not created directly
                OutputFile(output_template % "info"),
                OutputFile(output_template % "result"),
                OutputFile(output_template % "log"),
                # Only generated by newer versions of ExaML
                OutputFile(
                    os.path.basename(output_template % "modelFile"),
                    temporary=True,
                ),
            ],
            requirements=[EXAML_VERSION],
        )

        command.append_options(
            {
                # Ensures that output is saved to the temporary directory
                "-w": "%(TEMP_DIR)s",
                "-s": InputFile(input_binary),
                "-t": InputFile(initial_tree),
                "-n": "Pypeline",
                # Use the GAMMA model of NT substitution by default
                "-m": model,
            }
        )

        self._dirname = os.path.dirname(output_template)
        self._template = os.path.basename(output_template)

        CommandNode.__init__(
            self,
            command=command,
            description="running ExaML on %s" % (input_binary,),
            threads=threads,
            dependencies=dependencies,
        )

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
        checkpoints = glob.glob(os.path.join(temp, "ExaML_binaryCheckpoint.Pypeline_*"))
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
                    source = os.path.join(temp, filename)
                    destination = os.path.join(temp, self._template % match.groups())

                    fileutils.move_file(source, destination)

        CommandNode._teardown(self, config, temp)
