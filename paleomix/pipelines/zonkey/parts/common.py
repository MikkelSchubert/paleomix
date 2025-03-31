#!/usr/bin/python3
#
# Copyright (c) 2016 Mikkel Schubert <MikkelSch@gmail.com>
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
import paleomix.common.fileutils as fileutils
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
            description="writing sample-list to %s" % (output_file,),
            input_files=(config.database.filename,),
            output_files=(output_file,),
            dependencies=dependencies,
        )

    def _run(self, config, temp):
        (output_file,) = self.output_files
        samples = self._samples

        group = self._groups[max(self._groups)]
        group_colors = dict(zip(sorted(set(group.values())), _DEFAULT_COLORS))

        with open(fileutils.reroot_path(temp, output_file), "w") as handle:
            handle.write("Name\tGroup\tColor\n")

            for sample_name in sorted(samples):
                group_name = group[sample_name]
                group_color = group_colors[group_name]

                handle.write("%s\t%s\t%s\n" % (sample_name, group_name, group_color))

            handle.write("Sample\t-\t#000000\n")

    def _teardown(self, config, temp):
        (destination,) = self.output_files
        source = fileutils.reroot_path(temp, destination)

        fileutils.move_file(source, destination)
