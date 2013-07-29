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
from collections import defaultdict

from pypeline.common.sequences import split
from pypeline.common.fileutils import open_ro
from pypeline.common.formats.fasta import FASTA, FASTAError


class MSAError(FASTAError):
    pass


class MSA(frozenset):
    """Represents a Multiple Sequence Alignment of FASTA records."""

    def __new__(cls, sequences):
        records, names = [], set()
        for record in sequences:
            if record.name in names:
                raise MSAError("Duplicate name found in FASTA records: %r" % record.name)
            records.append(record)
            names.add(record.name)

        instance = frozenset.__new__(cls, records)
        MSA.validate(instance)
        return instance

    def seqlen(self):
        """Retrurns the length of the sequences in the MSA."""
        return len(iter(self).next().sequence)

    def exclude(self, names):
        """Builds a new MSA that excludes the named set of records."""
        included, excluded = [], set(names)
        for record in self:
            if record.name not in names:
                included.append(record)
            else:
                excluded.remove(record.name)
        if excluded:
            raise KeyError("Key(s) not found: %r" % (", ".join(map(str, excluded))))
        return MSA(included)


    def split(self, split_by = "123"):
        """Splits a MSA and returns a dictionary of keys to MSAs,
        using the keys in the 'split_by' parameter at the top
        level. See also pypeline.common.sequences.split."""
        self.validate(self)
        if not split_by:
            raise TypeError("No partitions to split by specified")

        results = dict((key, set()) for key in split_by)
        for record in self:
            for (key, partition) in split(record.sequence, split_by).iteritems():
                results[key].add(FASTA(record.name, None, partition))

        for (key, value) in results.items():
            results[key] = MSA(value)

        return results

    @classmethod
    def join(cls, *msas):
        """Merge multiple MSAs into a single MSA, by concatenating sequences in
        the order of the passed MSAs. Sequences are joined by name, and all MSAs
        must therefore contain the same set of sequence names. Meta information
        is not preserved."""
        cls.validate(*msas)

        merged = defaultdict(list)
        for msa in msas:
            for record in msa:
                merged[record.name].append(record.sequence)

        sequences = []
        for (name, sequence) in merged.iteritems():
            sequences.append(FASTA(name, None, "".join(sequence)))
        return MSA(sequences)


    @classmethod
    def from_lines(cls, lines):
        """Parses a MSA from a file/list of lines, and returns a dictionary
        of names to sequences. If read_meta is True, meta information included
        after the first space in header of each sequence:
          >NAME META-INFORMATION
          SEQUENCE
        As suggested above, sequences are expected to be in FASTA format."""
        return MSA(FASTA.from_lines(lines))

    @classmethod
    def from_file(cls, filename):
        """Reads a MSA from the specified filename. The file may
        be uncompressed, gzipped or bzipped. See also 'MSA.from_lines'."""
        fasta_file = open_ro(filename)
        try:
            return MSA.from_lines(fasta_file)
        finally:
            fasta_file.close()

    def to_file(self, fileobj):
        for fst in sorted(self):
            fileobj.write(str(fst))

    @classmethod
    def validate(cls, *msas):
        """Validates one or more MSAs, requiring:
        1. Identical sets of sequence names across all MSAs.
        2. That all names are non-empty strings.
        3. That all sequences are of the same length (per MSA).
        4. That no empty MSA (no sequences) are specified."""
        if not msas:
            raise TypeError("No MSAs given as arguments")

        seqs_all = msas[0].names()
        seqs_common = set(seqs_all)
        for msa in msas:
            if len(set(len(record.sequence) for record in msa)) != 1:
                raise MSAError("MSA contains sequences of differing lengths")

            seqs_all.update(msa.names())
            seqs_common &= set(msa.names())

        if seqs_all != seqs_common:
            raise MSAError("Some sequences not found in all MSAs: '%s'" \
                           % ("', '".join(seqs_all - seqs_common),))

    def __repr__(self):
        def _fasta_to_str(fst):
            return "FASTA(%r, %r, %r)" % \
              (fst.name, fst.meta, fst.sequence)
        return "MSA(%s)" % (", ".join(map(_fasta_to_str, sorted(self))))

    def names(self):
        return set(record.name for record in self)
