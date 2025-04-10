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

from paleomix.common import rtools
from paleomix.common.command import InputFile, OutputFile, TempOutputFile
from paleomix.common.fileutils import PathTypes
from paleomix.common.formats.newick import Newick
from paleomix.node import CommandNode
from paleomix.tools import factory


class MitoConsensusNode(CommandNode):
    def __init__(self, database, bamfile, output_prefix, dependencies=()):
        command = factory.new(
            [
                "zonkey:mito",
                InputFile(database),
                InputFile(bamfile),
                TempOutputFile(output_prefix),
            ],
            extra_files=[
                OutputFile(output_prefix + ".phy"),
                OutputFile(output_prefix + ".fasta"),
                OutputFile(output_prefix + ".summary"),
            ],
        )

        CommandNode.__init__(
            self,
            description=f"building consensus mitochondria from {bamfile}",
            command=command,
            dependencies=dependencies,
        )


class DrawPhylogenyNode(CommandNode):
    def __init__(self, samples, treefile, bootstraps, output_prefix, dependencies=()):
        command = factory.rscript(
            (
                os.path.join("zonkey", "tinytree.r"),
                # Temporary file generated in _setup
                TempOutputFile("rerooted.newick"),
                InputFile(samples),
                TempOutputFile(output_prefix),
            ),
            extra_files=[
                InputFile(bootstraps),
                OutputFile(output_prefix + ".pdf"),
                OutputFile(output_prefix + ".png"),
            ],
            requirements=[
                rtools.requirement("ape"),
                rtools.requirement("ggplot2"),
                rtools.requirement("grid"),
            ],
        )

        self._treefile = treefile
        self._bootstraps = bootstraps

        CommandNode.__init__(
            self,
            description=f"drawing phylogeny in {treefile}",
            command=command,
            dependencies=dependencies,
        )

    def _setup(self, temp: PathTypes) -> None:
        with open(self._bootstraps) as handle:
            bootstraps = [Newick.from_string(line.strip()) for line in handle]

        with open(self._treefile) as handle:
            tree = Newick.from_string(handle.read().strip())

        tree = tree.reroot_on_midpoint()
        tree = tree.add_support(bootstraps, "{Percentage:.0f}")
        with open(os.path.join(temp, "rerooted.newick"), "w") as handle:
            handle.write(f"{tree}\n")

        super()._setup(temp)
