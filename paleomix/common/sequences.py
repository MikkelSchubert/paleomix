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
from typing import Dict, Iterable, List


def _build_complement_table() -> str:
    # Pairs of complementary bases and ambiguous bases
    table = ["N"] * 256
    for nuc_a, nuc_b in ["AT", "CG", "NN", "RY", "KM", "SS", "WW", "BV", "DH", "XX"]:
        # Complement both upper/lower-case bases
        for _func in (str.upper, str.lower):
            table[ord(_func(nuc_a))] = _func(nuc_b)
            table[ord(_func(nuc_b))] = _func(nuc_a)

    return "".join(table)


IUPAC_TABLE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "N": "N",
    "AG": "R",
    "CT": "Y",
    "GT": "K",
    "AC": "M",
    "CG": "S",
    "AT": "W",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "N",
}


# Table for complementing sequences using `str.translate`
NT_COMPLEMENTS = _build_complement_table()
# IUPAC code to A, C, G, T nucleotides.
NT_CODES = dict(zip(IUPAC_TABLE.values(), IUPAC_TABLE))


def complement(sequence: str) -> str:
    """Returns the complement of a DNA sequence (string)."""
    return sequence.translate(NT_COMPLEMENTS)


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    return complement(sequence)[::-1]


def encode_genotype(nucleotides: Iterable[str]) -> str:
    """Parses a string representing a set of nucleotides observed at a loci,
    and returns the corresponding IUPAC code. Commas are allowed, but are
    simply ignored if found in the string. Does not handle lower-case
    nucleotides, due to lack of clear criteria for mixed case input.
    See e.g. http://www.ebi.ac.uk/2can/tutorials/aa.html"""
    chars = set(nucleotides).difference(",")
    if chars.difference("ACGTN"):
        raise ValueError(nucleotides)

    return IUPAC_TABLE["".join(sorted(chars))]


def split(sequence: str, split_by: str = "123") -> Dict[str, str]:
    """Splits a sequence by position, as specified by the 'split_by' parameter. By
    default, the function will split by codon position, and return a dictionary
    containing the keys '1', '2' and '3'.

    The 'split_by' parameter may contain any non-zero number of values, which must
    however be hashable. If a value is specified multiple times, then those positions
    are interleaved (e.g. split_by = "112" returns the first two positions in a codon
    as one sequence, as well as the last positions as one sequence."""
    if not split_by:
        raise ValueError("No split_by specified")

    results: Dict[str, List[str]] = {key: [] for key in split_by}
    keys = itertools.chain(itertools.cycle(split_by))
    for key, nucleotide in zip(keys, sequence):
        results[key].append(nucleotide)

    return {key: "".join(value) for key, value in results.items()}
