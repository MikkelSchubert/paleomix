#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import os
import random
import re
from collections.abc import Iterable

from paleomix.common import fileutils, versions
from paleomix.common.command import (
    AtomicCmd,
    InputFile,
    OptionsType,
    OutputFile,
    TempOutputFile,
)
from paleomix.node import CommandNode, Node

RAXML_VERSION = versions.Requirement(
    call=("raxmlHPC", "-version"),
    regexp=r"version (\d+\.\d+\.\d+)",
    specifiers=">=8.2.9",
)
RAXML_PTHREADS_VERSION = versions.Requirement(
    call=("raxmlHPC-PTHREADS", "-version"),
    regexp=r"version (\d+\.\d+\.\d+)",
    specifiers=">=8.2.9",
)


class RAxMLRapidBSNode(CommandNode):
    def __init__(
        self,
        *,
        input_alignment: str,
        output_template: str,
        input_partition: str | None = None,
        model: str = "GTRGAMMAI",
        replicates: str | int = "autoMRE",
        threads: int = 1,
        dependencies: Iterable[Node] = (),
        rng: random.Random | None = None,
    ) -> None:
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
        self._symlinks = [input_alignment, input_partition]
        self._template = os.path.basename(output_template)
        if rng is None:
            rng = random.Random()

        options: OptionsType = {
            # Perform rapid bootstrapping
            "-f": "a",
            # Output files are saved with a .PALEOMIX postfix, and subsequently renamed
            "-n": "PALEOMIX",
            # Ensures that output is saved to the temporary directory
            "-w": "%(TEMP_DIR)s",
            # Use the GTRGAMMA model of NT substitution by default
            "-m": model,
            # Enable Rapid Bootstrapping and set random seed
            "-x": rng.randint(0, 2**31 - 1),
            # Set random seed for parsimony inference
            "-p": rng.randint(0, 2**31 - 1),
            # Terminate bootstrapping upon convergence, not after N repetitions
            "-N": replicates,
        }

        extra_files = [
            # Input files, are not used directly (see below)
            InputFile(input_alignment),
            # Final output files, are not created directly
            OutputFile(output_template % "info"),
            OutputFile(output_template % "bestTree"),
            OutputFile(output_template % "bootstrap"),
            OutputFile(output_template % "bipartitions"),
            OutputFile(output_template % "bipartitionsBranchLabels"),
        ]

        # Symlink to sequence and partitions, to prevent the creation of *.reduced files
        # outside temp folder. In addition, it may be necessary to remove the .reduced
        # files if created
        alignment_basename = os.path.basename(input_alignment)
        options["-s"] = TempOutputFile(alignment_basename)
        extra_files.append(TempOutputFile(alignment_basename + ".reduced"))

        if input_partition is not None:
            partition_basename = os.path.basename(input_partition)

            options["-q"] = TempOutputFile(partition_basename)
            extra_files.append(InputFile(input_partition))
            extra_files.append(TempOutputFile(partition_basename + ".reduced"))

        if threads > 1:
            command = AtomicCmd(
                ["raxmlHPC-PTHREADS", "-T", threads],
                extra_files=extra_files,
                requirements=[RAXML_PTHREADS_VERSION],
            )
        else:
            command = AtomicCmd(
                "raxmlHPC",
                extra_files=extra_files,
                requirements=[RAXML_VERSION],
            )

        command.append_options(options)

        CommandNode.__init__(
            self,
            command=command,
            description=f"inferring phylogeny from {input_alignment} using RAxML",
            threads=threads,
            dependencies=dependencies,
        )

    def _setup(self, temp: fileutils.PathTypes) -> None:
        super()._setup(temp)

        # Required to avoid the creation of files outside the temp folder
        for filename in self._symlinks:
            if filename is not None:
                source = os.path.abspath(filename)
                destination = os.path.join(temp, os.path.basename(filename))

                os.symlink(source, destination)

    def _run(self, temp: fileutils.PathTypes) -> None:
        # RAxML needs to be run with an absolute path for -w
        super()._run(os.path.abspath(temp))

    def _teardown(self, temp: fileutils.PathTypes) -> None:
        for filename in os.listdir(temp):
            match = re.match("RAxML_(.*).PALEOMIX", filename)
            if match:
                source = os.path.join(temp, filename)
                destination = os.path.join(temp, self._template % match.groups())

                fileutils.move_file(source, destination)

        super()._teardown(temp)
