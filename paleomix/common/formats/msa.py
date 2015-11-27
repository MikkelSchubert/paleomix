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
from itertools import izip
from collections import defaultdict

from paleomix.common.sequences import split
from paleomix.common.fileutils import open_ro
from paleomix.common.formats.fasta import FASTA, FASTAError
from paleomix.common.sequences import NT_CODES, encode_genotype
from paleomix.common.utilities import safe_coerce_to_frozenset


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

        if not records:
            raise MSAError("MSA does not contain any sequences")

        instance = frozenset.__new__(cls, records)
        MSA.validate(instance)
        return instance

    def seqlen(self):
        """Retrurns the length of the sequences in the MSA."""
        return len(iter(self).next().sequence)

    def exclude(self, names):
        """Builds a new MSA that excludes the named set of records."""
        _, excluded, _ = self._group(names)
        return MSA(excluded)

    def select(self, names):
        """Builds a new MSA that includes only the named set of records."""
        included, _, _ = self._group(names)
        return MSA(included)

    def reduce(self):
        columns = []
        uncalled = frozenset("Nn-")
        for column in izip(*(record.sequence for record in self)):
            if (frozenset(column) - uncalled):
                columns.append(column)

        if not columns:
            return None

        records = []
        for (record, sequence) in izip(self, izip(*columns)):
            records.append(FASTA(record.name, record.meta, "".join(sequence)))

        return MSA(records)

    def filter_singletons(self, to_filter, filter_using):
        included, excluded, to_filter \
            = self._group(filter_using, to_filter)

        sequence = list(to_filter.sequence)
        sequences = [record.sequence.upper() for record in included]
        for (index, nts) in enumerate(zip(*sequences)):
            current_nt = sequence[index].upper()
            if current_nt in "N-":
                continue

            allowed_nts = set()
            for allowed_nt in nts:
                if allowed_nt not in "N-":
                    allowed_nts.update(NT_CODES[allowed_nt])
            filtered_nts = frozenset(NT_CODES[current_nt]) & allowed_nts

            if not filtered_nts:
                filtered_nts = "N"

            genotype = encode_genotype(filtered_nts)
            if genotype != current_nt:
                sequence[index] = genotype.lower()
        new_record = FASTA(to_filter.name,
                           to_filter.meta,
                           "".join(sequence))

        return MSA([new_record] + included + excluded)


    def split(self, split_by = "123"):
        """Splits a MSA and returns a dictionary of keys to MSAs,
        using the keys in the 'split_by' parameter at the top
        level. See also paleomix.common.sequences.split."""
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
        except MSAError, error:
            raise MSAError("%s in file %r" % (error, filename))
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


    def _group(self, selection, extra = None):
        selection = safe_coerce_to_frozenset(selection)
        if (extra in selection):
            raise MSAError("Key used for multiple selections: %r" % extra)
        elif not selection:
            raise ValueError("No FASTA names given")

        missing_keys = selection - self.names()
        if missing_keys:
            raise KeyError("Key(s) not found: %r" % (", ".join(map(str, missing_keys))))

        included, excluded, other = [], [], None
        for record in self:
            if record.name in selection:
                included.append(record)
            elif record.name != extra:
                excluded.append(record)
            else:
                other = record

        return included, excluded, other
