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
"""Various functions relating to DNA sequence manipulation."""

import itertools


# Pairs of complementary bases and ambigious basees
_COMPL = ["AT", "CG", "NN", "RY", "KM", "SS", "WW", "BV", "DH", "XX"]
_COMPL_TABLE = ["N"] * 256
for (_a, _b) in _COMPL:
    # Complement both upper/lower-case bases
    for _func in (str.upper, str.lower):
        _COMPL_TABLE[ord(_func(_a))] = _func(_b)
        _COMPL_TABLE[ord(_func(_b))] = _func(_a)
_COMPL_TABLE = "".join(_COMPL_TABLE)


# Table of nt codes (IUPAC codes) used to encode (ambigious) bases:
#   Nomenclature for incompletely specified bases in nucleic acid sequences.
#   Recommendations 1984. J Biol Chem. 1986 Jan 5;261(1):13-7.
#   PubMed PMID: 2416744.
NT_CODES = [
    ["A", "A"],
    ["C", "C"],
    ["G", "G"],
    ["T", "T"],
    ["N", "N"],
    ["R", "AG"],
    ["Y", "CT"],
    ["K", "GT"],
    ["M", "AC"],
    ["S", "CG"],
    ["W", "AT"],
    ["B", "CGT"],
    ["D", "AGT"],
    ["H", "ACT"],
    ["V", "ACG"],
    ["N", "ACGT"],
]

_NT_CODES_TABLE = {}
for (_abr, _nts) in NT_CODES:
    _NT_CODES_TABLE[frozenset(_nts)] = _abr
    _NT_CODES_TABLE[frozenset(_nts + ",")] = _abr

NT_CODES = dict(NT_CODES)


def complement(sequence):
    """Returns the complement of a DNA sequence (string)."""
    return sequence.translate(_COMPL_TABLE)


def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence."""
    return complement(sequence)[::-1]


def encode_genotype(nucleotides):
    """Parses a string representing a set of nucleotides observed at a loci,
    and returns the corresponding IUPAC code. Commas are allowed, but are
    simply ignored if found in the string. Does not handle lower-case
    nucleotides, due to lack of clear criteria for mixed case input.
    See e.g. http://www.ebi.ac.uk/2can/tutorials/aa.html"""
    try:
        return _NT_CODES_TABLE[frozenset(nucleotides)]
    except KeyError:
        raise ValueError(nucleotides)


def split(sequence, split_by="123"):
    """Splits a sequence by position, as specified by the 'split_by' parameter. By
    default, the function will split by codon position, and return a dictionary
    containing the keys '1', '2' and '3'.

    The 'split_by' parameter may contain any non-zero number of values, which must
    however be hashable. If a value is specified multiple times, then those positions
    are interleaved (e.g. split_by = "112" returns the first two positions in a codon
    as one sequence, as well as the last positions as one sequence."""
    if not split_by:
        raise ValueError("No split_by specified")

    results = dict((key, []) for key in split_by)
    keys = itertools.chain(itertools.cycle(split_by))
    for (key, nucleotide) in zip(keys, sequence):
        results[key].append(nucleotide)

    for key in results:
        results[key] = "".join(results[key])

    return results
