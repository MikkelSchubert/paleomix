# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

from paleomix.common import fileutils
from paleomix.node import Node

_DEFAULT_COLORS = (
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
)


class WriteSampleList(Node):
    def __init__(self, config, output_file, dependencies=()):
        self._samples = config.database.samples
        self._groups = config.database.groups

        Node.__init__(
            self,
            description=f"writing sample-list to {output_file}",
            input_files=(config.database.filename,),
            output_files=(output_file,),
            dependencies=dependencies,
        )

    def _run(self, temp):
        (output_file,) = self.output_files
        samples = self._samples

        group = self._groups[max(self._groups)]
        group_colors = dict(zip(sorted(set(group.values())), _DEFAULT_COLORS))

        with open(fileutils.reroot_path(temp, output_file), "w") as handle:
            handle.write("Name\tGroup\tColor\n")

            for sample_name in sorted(samples):
                group_name = group[sample_name]
                group_color = group_colors[group_name]

                handle.write(f"{sample_name}\t{group_name}\t{group_color}\n")

            handle.write("Sample\t-\t#000000\n")

    def _teardown(self, temp: fileutils.PathTypes) -> None:
        (destination,) = self.output_files
        source = fileutils.reroot_path(temp, destination)

        fileutils.move_file(source, destination)
