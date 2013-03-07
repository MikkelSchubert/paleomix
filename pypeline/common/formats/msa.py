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
import gzip
import types
from collections import defaultdict

from pypeline.common.sequences import split
from pypeline.common.formats.fasta import parse_fasta, print_fasta


class MSAError(RuntimeError):
    pass


def split_msa(msa, split_by = "123"):
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
    validate_msa(*msas)
    results = defaultdict(list)
    for msa in msas:
        for (name, sequence) in msa.iteritems():
            results[name].append(sequence)

    return dict((key, "".join(value)) for (key, value) in results.iteritems())
    

def parse_msa(lines, read_header = False):
    msa, headers = {}, {}
    for (header, sequence) in parse_fasta(lines):
        if not header:
            raise MSAError("MSA record without name found")

        name = header.split(None, 1)[0]
        if name in msa:
            raise MSAError("Duplicate names found, cannot be represented as MSA")
        msa[name] = sequence
        headers[name] = header

    validate_msa(msa)
    if read_header:
        return msa, header
    return msa


def read_msa(filename, read_header = False):
    func = gzip.open if filename.endswith(".gz") else open
    fasta_file = func(filename, "r")
    try:
        return parse_msa(iter(fasta_file), read_header = read_header)
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
    if not msas:
        raise TypeError("No MSAs given as arguments")

    keywords = set(msas[0])
    for msa in msas:
        if not msa:
            raise MSAError("MSA with no sequences found")
        elif not all((name and isinstance(name, types.StringTypes)) for name in msa):
            raise MSAError("Names in MSA must be non-empty strings")
        elif len(set(len(seq) for seq in msa.itervalues())) != 1:
            raise MSAError("MSA contains sequences of differing lengths")
        elif set(msa) != keywords:
            raise MSAError("MSAs contain mismatching sequences")
