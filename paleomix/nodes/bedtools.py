#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
from paleomix.node import Node, NodeError

from paleomix.common.bedtools import read_bed_file, pad_bed_records, merge_bed_records
from paleomix.common.fileutils import move_file, reroot_path


class PaddedBedNode(Node):
    """Simple node for padding BED records a fixed amount and merging
    overlapping records. Columns beyond the 3rd column are dropped.
    """

    def __init__(self, infile, outfile, fai_file, amount=0, dependencies=()):
        self._amount = int(amount)
        self._infile = infile
        self._outfile = outfile
        self._fai_file = fai_file

        Node.__init__(
            self,
            description="padding BED records (%i) in %s" % (amount, infile),
            input_files=(infile, fai_file),
            output_files=(outfile,),
            dependencies=dependencies,
        )

    def _run(self, config, temp):
        contigs = {}
        with open(self._fai_file) as handle:
            for line in handle:
                name, length, _ = line.split("\t", 2)
                if name in contigs:
                    raise NodeError(
                        "Reference genome contains multiple "
                        "identically named contigs (%r)!" % (name,)
                    )

                contigs[name] = int(length)

        with open(reroot_path(temp, self._outfile), "w") as handle:
            records = list(read_bed_file(self._infile, contigs=contigs))
            pad_bed_records(records=records, padding=self._amount, max_sizes=contigs)

            for record in merge_bed_records(records):
                handle.write("%s\n" % (record,))

    def _teardown(self, config, temp):
        source = reroot_path(temp, self._outfile)
        move_file(source, self._outfile)
