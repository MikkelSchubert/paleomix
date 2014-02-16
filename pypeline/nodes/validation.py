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
import collections
from itertools import \
    izip_longest

import pysam

from pypeline.node import \
    Node, \
    NodeError
from pypeline.common.fileutils import \
    describe_files
from pypeline.common.utilities import \
    chain_sorted


class DetectInputDuplicationNode(Node):
    """Attempts to detect reads included multiple times as input based on the
    presence of reads with identical names AND sequences. This is compromise
    between sensitivity, specificity, and running time.

    A possible refinement would be to consider reads with the same name where
    one read is the prefix of the other (due to different amounts of trimming
    or collapsing of reads).
    """

    def __init__(self, input_files, output_file, dependencies=()):
        Node.__init__(self,
                      description="<DetectInputDuplication: %s>"
                      % (describe_files(input_files)),
                      input_files=input_files,
                      output_files=output_file,
                      dependencies=dependencies)

    def run(self, _):
        handles = []
        try:
            sequences = []
            for fpath in self.input_files:
                handle = pysam.Samfile(fpath)
                handles.append(handle)

                sequence = izip_longest(handle, (), fillvalue=fpath)
                sequences.append(sequence)

            position = 0
            records = chain_sorted(*sequences, key=self._key_by_tid_pos)
            observed_reads = collections.defaultdict(list)
            for (record, fpath) in records:
                if record.pos != position:
                    self._process_reads(observed_reads)
                    observed_reads.clear()
                    position = record.pos

                key = (record.qname, record.seq)
                observed_reads[key].append(fpath)
            self._process_reads(observed_reads)

            # Everything is ok, touch the output files
            for fpath in self.output_files:
                with open(fpath, "w"):
                    pass
        finally:
            for handle in handles:
                handle.close()

    @classmethod
    def _process_reads(cls, observed_reads):
        for ((name, _), fpaths) in observed_reads.iteritems():
            if len(fpaths) > 1:
                message = ["Read %r found in multiple files:" % (name,)]
                for fpath in fpaths:
                    message.append("  - %r" % (fpath,))
                message.append("")
                message.append("This indicates that the same data files have "
                               "been included multiple times in the project. "
                               "Please review the input files used in this "
                               "project, to ensure that each set of data is "
                               "included only once.")
                raise NodeError("\n".join(message))

    @classmethod
    def _key_by_tid_pos(cls, record):
        return (record[0].tid, record[0].pos)
