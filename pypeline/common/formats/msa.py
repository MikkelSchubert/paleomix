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
"""Functions for manipulating Multiple Sequence Alignments."""
import sys
import types
from collections import defaultdict

from pypeline.common.sequences import split
from pypeline.common.fileutils import open_ro
from pypeline.common.formats.fasta import parse_fasta, print_fasta, FASTAError


class MSAError(FASTAError):
    pass


def split_msa(msa, split_by = "123"):
    """Splits a MSA and returns a dictionary of keys to MSAs,
    using the keys in the 'split_by' parameter at the top
    level. See also pypeline.common.sequences.split."""
    validate_msa(msa)
    if not split_by:
        raise TypeError("No partitions to split by specified")

    results = {}
    for key in split_by:
        results[key] = dict((name, {}) for name in msa)

    for (name, sequence) in msa.iteritems():
        for (key, partition) in split(sequence, split_by).iteritems():
            results[key][name] = partition

    return results


def join_msa(*msas):
    """Merge multiple MSAs into a single MSA, by concatenating sequences in
    the order of the passed MSAs. Sequences are joined by name, and all MSAs
    must therefore contain the same set of sequence names."""
    validate_msa(*msas)
    results = defaultdict(list)
    for msa in msas:
        for (name, sequence) in msa.iteritems():
            results[name].append(sequence)

    return dict((key, "".join(value)) for (key, value) in results.iteritems())


def parse_msa(lines, read_meta = False):
    """Parses a MSA from a file/list of lines, and returns a dictionary
    of names to sequences. If read_meta is True, meta information included
    after the first space in header of each sequence:
      >NAME META-INFORMATION
      SEQUENCE
    As suggested above, sequences are expected to be in FASTA format."""
    msa, metas = {}, {}
    for ((name, meta), sequence) in parse_fasta(lines):
        if name in msa:
            raise MSAError("Duplicate names found, cannot be represented as MSA: " + name)
        msa[name] = sequence
        metas[name] = meta

    validate_msa(msa)
    if read_meta:
        return msa, metas
    return msa


def read_msa(filename, read_meta = False):
    """Reads a MSA from the specified filename. The file may
    be uncompressed, gzipped or bzipped. See also 'parse_msa'."""
    fasta_file = open_ro(filename)
    try:
        return parse_msa(iter(fasta_file), read_meta = read_meta)
    finally:
        fasta_file.close()


def print_msa(msa, file = sys.stdout):
    validate_msa(msa)
    for group in sorted(msa):
        print_fasta(group, msa[group], file)


def write_msa(msa, filename):
    validate_msa(msa)
    with open(filename, "w") as fileobj:
        print_msa(msa, fileobj)


def validate_msa(*msas):
    """Validates one or more MSAs, requiring:
     1. Identical sets of sequence names across all MSAs.
     2. That all names are non-empty strings.
     3. That all sequences are of the same length (per MSA).
     4. That no empty MSA (no sequences) are specified."""
    if not msas:
        raise TypeError("No MSAs given as arguments")

    seqs_all = set(msas[0])
    seqs_common = set(msas[0])
    for msa in msas:
        if not msa:
            raise MSAError("MSA with no sequences found")
        elif not all((name and isinstance(name, types.StringTypes)) for name in msa):
            raise MSAError("Names in MSA must be non-empty strings")
        elif len(set(len(seq) for seq in msa.itervalues())) != 1:
            raise MSAError("MSA contains sequences of differing lengths")

        seqs_all.update(msa)
        seqs_common &= set(msa)

    if seqs_all != seqs_common:
        raise MSAError("Some sequences not found in all MSAs: '%s'" \
                       % ("', '".join(seqs_all - seqs_common),))
