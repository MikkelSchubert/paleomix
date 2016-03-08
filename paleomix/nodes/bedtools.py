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
import paleomix.common.versions as versions

from paleomix.atomiccmd.command import \
    AtomicCmd
from paleomix.node import \
    CommandNode, \
    Node, \
    NodeError

from paleomix.common.bedtools import \
    read_bed_file
from paleomix.common.fileutils import \
    move_file, \
    reroot_path


BEDTOOLS_VERSION \
    = versions.Requirement(call=("bedtools", "--version"),
                           search=r"bedtools v?(\d+)\.(\d+)\.(\d+)",
                           checks=versions.GE(2, 15, 0))


class SlopBedNode(CommandNode):
    def __init__(self, infile, outfile, genome, from_start=0, from_end=0,
                 strand_relative=False, dependencies=()):
        if type(from_start) != type(from_end):
            raise ValueError("Parameters 'from_start' and 'from_end' should "
                             "be of same type!")

        call = ["bedtools", "slop",
                "-i", "%(IN_FILE)s",
                "-g", "%(IN_GENOME)s",
                "-l", str(from_start),
                "-r", str(from_end)]

        if strand_relative:
            call.append("-s")
        if type(from_start) is float:
            call.append("-pct")

        command = AtomicCmd(call,
                            IN_FILE=infile,
                            IN_GENOME=genome,
                            OUT_STDOUT=outfile,
                            CHECK_VERSION=BEDTOOLS_VERSION)

        description = "<SlopBed: '%s' -> '%s'>" % (infile, outfile)

        CommandNode.__init__(self,
                             description=description,
                             command=command,
                             dependencies=dependencies)


class PaddedBedNode(Node):
    """Simple node for padding BED records a fixed amount and merging
    overlapping records. Columns beyond the 3rd column are dropped.
    """

    def __init__(self, infile, outfile, fai_file, amount=0, dependencies=()):
        self._amount = int(amount)
        self._infile = infile
        self._outfile = outfile
        self._fai_file = fai_file

        Node.__init__(self,
                      description='<PaddedBed (%i): %r -> %r>'
                      % (amount, infile, outfile),
                      input_files=(infile, fai_file),
                      output_files=(outfile,),
                      dependencies=dependencies)

    def _run(self, config, temp):
        contigs = {}
        with open(self._fai_file) as handle:
            for line in handle:
                name, length, _ = line.split('\t', 2)
                if name in contigs:
                    raise NodeError('Reference genome contains multiple '
                                    'identically named contigs (%r)!'
                                    % (name,))

                contigs[name] = int(length)

        with open(reroot_path(temp, self._outfile), 'w') as handle:
            for record in read_bed_file(self._infile, contigs=contigs):
                max_length = contigs[record.contig]
                record.start = max(0, record.start - self._amount)
                record.end = min(record.end + self._amount, max_length)

                handle.write('%s\n' % (record,))

    def _teardown(self, config, temp):
        source = reroot_path(temp, self._outfile)
        move_file(source, self._outfile)
