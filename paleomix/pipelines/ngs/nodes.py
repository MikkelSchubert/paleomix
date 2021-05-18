#!/usr/bin/env python3
import os

import paleomix.common.rtools as rtools

from paleomix.node import CommandNode
from paleomix.atomiccmd.command import AtomicCmd, InputFile, AuxilleryFile, OutputFile


class TranchesPlotsNode(CommandNode):
    def __init__(self, input_table, output_prefix, dependencies=()):
        command = AtomicCmd(
            (
                "Rscript",
                AuxilleryFile(rtools.rscript("ngs", "tranches.r")),
                InputFile(input_table),
                OutputFile(os.path.basename(output_prefix), temporary=True),
            ),
            extra_files=(
                OutputFile(output_prefix + ".pdf"),
                OutputFile(output_prefix + ".txt"),
            ),
        )

        CommandNode.__init__(
            self,
            command=command,
            description="plotting recalibration tranches for {!r}".format(input_table),
            dependencies=dependencies,
        )
