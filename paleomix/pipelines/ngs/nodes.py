#!/usr/bin/env python3
import paleomix.common.rtools as rtools
from paleomix.common.command import (
    AtomicCmd,
    AuxilleryFile,
    InputFile,
    OutputFile,
    TempOutputFile,
)
from paleomix.node import CommandNode


class TranchesPlotsNode(CommandNode):
    def __init__(self, input_table, output_prefix, dependencies=()):
        command = AtomicCmd(
            (
                "Rscript",
                AuxilleryFile(rtools.rscript("ngs", "tranches.r")),
                InputFile(input_table),
                TempOutputFile(output_prefix),
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
