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
import os
import sys
import types

import pysam

from paleomix.common.utilities import \
     fragment, \
     split_before, \
     Immutable, \
     TotallyOrdered
from paleomix.common.fileutils import open_ro
from paleomix.common.formats._common import FormatError


class FASTAError(FormatError):
    pass


class FASTA(TotallyOrdered, Immutable):
    def __init__(self, name, meta, sequence):
        if not (name and isinstance(name, types.StringTypes)):
            raise FASTAError("FASTA name must be a non-empty string")
        elif not (isinstance(meta, types.StringTypes) or (meta is None)):
            raise FASTAError("FASTA meta must be a string, or None")
        elif not isinstance(sequence, types.StringTypes):
            raise FASTAError("FASTA sequence must be a string")

        Immutable.__init__(self,
                           name=name,
                           meta=meta,
                           sequence=sequence)

    def write(self, fileobj=sys.stdout):
        """Prints a FASTA sequence (iterable), wrapping long sequences at 60
        characters."""
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
                raise FASTAError("FASTA record does not contain sequence: %s"
                                 % (name[1:],))

            # Split out any meta information
            name_and_meta = name[1:].split(None, 1)
            if len(name_and_meta) < 2:
                name_and_meta.append(None)
            name, meta = name_and_meta

            yield FASTA(name=name,
                        meta=meta,
                        sequence="".join(record[1:]))

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

    @classmethod
    def index_and_collect_contigs(cls, filename):
        """Creates an index (.fai; if it does not already exist) for a FASTA
        file using 'pysam', and returns a dictionary of {contig: length} listed
        in that file; if the .fai file can not be created, or if the FASTA file
        contains sequences with identical names, then a FASTAError is raised.
        """
        fai_filename = filename + ".fai"
        if not os.path.exists(fai_filename):
            if not os.access(os.path.dirname(filename), os.W_OK):
                message = \
                    "FASTA index is missing, but folder is\n" \
                    "not writable, so it cannot be created:\n" \
                    "  Filename = %s\n\n" \
                    "Either change permissions on the folder, or move\n" \
                    "the FASTA file to different location." % (filename,)
                raise FASTAError(message)

            # Use pysam to index the file
            pysam.Fastafile(filename).close()

        names = set()
        contigs = []
        with open(fai_filename) as faihandle:
            for line in faihandle:
                name, length, _ = line.split(None, 2)
                if name in names:
                    raise FASTAError("Reference contains multiple identically "
                                     "named sequences:\n  Path = %r\n  Name = "
                                     "%r\nPlease ensure that sequences have "
                                     "unique names" % (filename, name))
                names.add(name)
                contigs.append((name, int(length)))

        return contigs

    def __lt__(self, other):
        if not isinstance(other, FASTA):
            return NotImplemented

        return (self.name, self.meta, self.sequence) \
            < (other.name, other.meta, other.sequence)

    def __hash__(self):
        return hash((self.name, self.meta, self.sequence))

    def __repr__(self):
        """Returns string representation of FASTA sequence, using the standard,
        FASTA file format, wrapping long sequences at 60 characters.
        """
        name = self.name
        if self.meta:
            name = "%s %s" % (name, self.meta)
        return ">%s\n%s\n" % (name, "\n".join(fragment(60, self.sequence)))
