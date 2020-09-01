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
import sys

import pysam

from paleomix.common.utilities import fragment, split_before, Immutable, TotallyOrdered
from paleomix.common.fileutils import open_ro, fspath
from paleomix.common.formats._common import FormatError


class FASTAError(FormatError):
    pass


class FASTA(TotallyOrdered, Immutable):
    def __init__(self, name, meta, sequence):
        if not (name and isinstance(name, str)):
            raise FASTAError("FASTA name must be a non-empty string")
        elif not (isinstance(meta, str) or (meta is None)):
            raise FASTAError("FASTA meta must be a string, or None")
        elif not isinstance(sequence, str):
            raise FASTAError("FASTA sequence must be a string")

        Immutable.__init__(self, name=name, meta=meta or "", sequence=sequence)

    def write(self, fileobj=sys.stdout):
        """Writes a FASTA record to fileobj, wrapping sequences at 60 chars"""
        name = self.name
        if self.meta:
            name = "%s %s" % (name, self.meta)

        fileobj.write(">%s\n%s\n" % (name, "\n".join(fragment(60, self.sequence))))

    @classmethod
    def from_lines(cls, lines):
        """Parses FASTA sequences found in a sequence of lines, and returns
        a tuple for each FASTA record: ((name, meta-information), sequence)
        No assumptions are made about the line-lengths."""
        lines = (line.rstrip() for line in lines)
        for record in split_before(lines, lambda v: v.startswith(">")):
            name = record[0]
            if (not name.startswith(">")) or (len(name) == 1):
                raise FASTAError("Unnamed FASTA record")
            elif len(record) == 1:
                raise FASTAError(
                    "FASTA record does not contain sequence: %s" % (name[1:],)
                )

            # Split out any meta information
            name_and_meta = name[1:].split(None, 1)
            if len(name_and_meta) < 2:
                name_and_meta.append("")
            name, meta = name_and_meta

            yield FASTA(name=name, meta=meta, sequence="".join(record[1:]))

    @classmethod
    def from_file(cls, filename):
        """Reads an unindexed FASTA file, returning a sequence of
        tuples containing the name and sequence of each entry in
        the file. The FASTA file may be GZIP/BZ2 compressed."""
        with open_ro(filename) as fasta_file:
            yield from FASTA.from_lines(fasta_file)

    @classmethod
    def index_and_collect_contigs(cls, filename):
        """Creates an index (.fai; if it does not already exist) for a FASTA file using
        pysam, and returns a dictionary of {contig: length} listed in that file.
        """
        # If an index does not already exist, then it is created automatically
        with pysam.FastaFile(fspath(filename)) as handle:
            return dict(zip(handle.references, handle.lengths))

    def __lt__(self, other):
        if not isinstance(other, FASTA):
            return NotImplemented

        return (self.name, self.meta, self.sequence) < (
            other.name,
            other.meta,
            other.sequence,
        )

    def __hash__(self):
        return hash((self.name, self.meta, self.sequence))

    def __repr__(self):
        return "FASTA(%r, %r, %r)" % (self.name, self.meta, self.sequence)
