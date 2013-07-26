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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
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
import types

from pypeline.common.utilities import fragment, split_before
from pypeline.common.fileutils import open_ro
from pypeline.common.formats._common import FormatError


class FASTAError(FormatError):
    pass


class FASTA:
    def __init__(self, name, meta, sequence):
        self._name     = name
        self._meta     = meta
        self._sequence = sequence

        if not (name and isinstance(name, types.StringTypes)):
            raise FASTAError("FASTA name must be a non-empty string")
        elif not (isinstance(meta, types.StringTypes) or (meta is None)):
            raise FASTAError("FASTA meta-information must be a string, or None")
        elif not (sequence or isinstance(name, types.StringTypes)):
            raise FASTAError("FASTA sequence must be a non-empty string")

    @property
    def name(self):
        return self._name

    @property
    def meta(self):
        return self._meta

    @property
    def sequence(self):
        return self._sequence


    def write(self, fileobj = sys.stdout):
        """Prints a FASTA sequence (iterable), wrapping long sequences at 60 chars."""
        fileobj.write(repr(self))

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
                raise FASTAError("FASTA record does not contain sequence: " + name[1:])

            # Split out any meta information
            name_and_meta = name[1:].split(None, 1)
            if len(name_and_meta) < 2:
                name_and_meta.append(None)
            name, meta = name_and_meta

            yield FASTA(name     = name,
                        meta     = meta,
                        sequence = "".join(record[1:]))

    @classmethod
    def from_file(cls, filename):
        """Reads an unindexed FASTA file, returning a sequence of
        tuples containing the name and sequence of each entry in
        the file. The FASTA file may be GZIP/BZ2 compressed."""
        fasta_file = open_ro(filename)
        try:
            for record in FASTA.from_lines(fasta_file):
                yield record
        finally:
            fasta_file.close()

    def __repr__(self):
        """Process a printable FASTA sequence, wrapping long sequences at 60 chars."""
        name = self.name
        if self.meta:
            name = "%s %s" % (name, self.meta)
        return ">%s\n%s\n" % (name, "\n".join(fragment(60, self.sequence)))


    def __eq__(self, other):
        return (self.name, self.meta, self.sequence) == \
          (other.name, other.meta, other.sequence)

    def __ne__(self, other):
        return not (self == other)


    def __hash__(self):
        return hash((self._name, self._meta, self._sequence))
