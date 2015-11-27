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
"""Various functions relating to DNA sequence manipulation."""

import itertools


# Pairs of complementary bases and ambigious basees
_COMPL = [ "AT", "CG",
           "NN", "RY",
           "KM", "SS",
           "WW", "BV",
           "DH", "XX" ]
_COMPL_TABLE = ["N"] * 256
for (_a, _b) in _COMPL :
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
    ["N", "ACGT"]]

_NT_CODES_TABLE = {}
for (_abr, _nts) in NT_CODES:
    _NT_CODES_TABLE[frozenset(_nts)] = _abr
    _NT_CODES_TABLE[frozenset(_nts + ",")] = _abr

NT_CODES = dict(NT_CODES)


CODONS = {
    "+": {
        "TTT": "Phe",  "TCT": "Ser",  "TAT": "Tyr",  "TGT": "Cys",
        "TTC": "Phe",  "TCC": "Ser",  "TAC": "Tyr",  "TGC": "Cys",
        "TTA": "Leu",  "TCA": "Ser",  "TAA": "Stop", "TGA": "Stop",
        "TTG": "Leu",  "TCG": "Ser",  "TAG": "Stop", "TGG": "Trp",

        "CTT": "Leu",  "CCT": "Pro",  "CAT": "His",  "CGT": "Arg",
        "CTC": "Leu",  "CCC": "Pro",  "CAC": "His",  "CGC": "Arg",
        "CTA": "Leu",  "CCA": "Pro",  "CAA": "Gln",  "CGA": "Arg",
        "CTG": "Leu",  "CCG": "Pro",  "CAG": "Gln",  "CGG": "Arg",

        "ATT": "Ile",  "ACT": "Thr",  "AAT": "Asn",  "AGT": "Ser",
        "ATC": "Ile",  "ACC": "Thr",  "AAC": "Asn",  "AGC": "Ser",
        "ATA": "Ile",  "ACA": "Thr",  "AAA": "Lys",  "AGA": "Arg",
        "ATG": "Met",  "ACG": "Thr",  "AAG": "Lys",  "AGG": "Arg",

        "GTT": "Val",  "GCT": "Ala",  "GAT": "Asp",  "GGT": "Gly",
        "GTC": "Val",  "GCC": "Ala",  "GAC": "Asp",  "GGC": "Gly",
        "GTA": "Val",  "GCA": "Ala",  "GAA": "Glu",  "GGA": "Gly",
        "GTG": "Val",  "GCG": "Ala",  "GAG": "Glu",  "GGG": "Gly"},
    "-": {}}


def complement(sequence):
    """Returns the complement of a DNA sequence (string)."""
    return sequence.translate(_COMPL_TABLE)


def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence."""
    return complement(sequence)[::-1]


for _codon, _aa in CODONS["+"].iteritems():
    assert not set(_codon) - set("ACGT")
    CODONS["-"][reverse_complement(_codon)] = _aa

assert len(CODONS["+"]) == 64
assert len(CODONS["-"]) == 64


def encode_genotype(nucleotides):
    """Parses a string representing a set of nucleotides observed at a loci,
    and returns the corresponding IUPAC code. Commas are allowed, but are
    simply ignored if found in the string. Does not handle lower-case
    nucleotides, due to lack of clear criteria for mixed case input.
    See e.g. http://www.ebi.ac.uk/2can/tutorials/aa.html"""
    try:
        return _NT_CODES_TABLE[frozenset(nucleotides)]
    except KeyError:
        raise ValueError("Invalid input for 'encode_genotype': %s" % (repr(nucleotides), ))



def count_nts(sequence):
    """Given a nucleotide sequence (str), this function returns
    the number of each type of nucleotide representable using
    IUPAC codes. The sequence must not contain non-IUPAC
    nucleotides, or other annotation. IUPAC nucleotides are
    handled in a case-insensitive manner."""
    counts = {}
    sequence = sequence.upper()
    for nucleotide in NT_CODES:
        count = sequence.count(nucleotide)
        if count:
            counts[nucleotide] = count

    if len(sequence) != sum(counts.itervalues()):
        raise ValueError("Sequence contains non-(IUPAC-)nucleotides: %s" % \
                         ", ".join(set(sequence) - set(counts)))

    return counts


def count_gc_diploid(sequence):
    """Given a sequence, this function returns a tuple containing the
    the total number of bases that were G/C, as well as the total number
    of bases. The sequence is assumed to represent a diploid genome, with
    the total number of bases being twice the sequence length, and
    hence IUPAC codes representing one of or both of G/C are treated as
    reflecting both strands. Thus R counts for 1, while S counts for 2.

    The sequence must only contain valid IUPAC codes, and no other
    form of annotation. Both uppercase/lowercase G/Cs are counted.
    Ambigious site (n/N) are not counted, neither in the number of G/C,
    nor in the total number of bases."""
    total_nts = total_gc = 0
    counts = count_nts(sequence)
    for (code, count) in counts.iteritems():
        value = 0
        if code == "N":
            continue
        elif code in "CGcg":
            value = 2
        else:
            code_represents = NT_CODES[code]
            if (len(code_represents) > 2) and (code != 'N'):
                raise ValueError("calculate_gcp assumes diploid genome, nt code for tri-valued SNP observed: " + code)

            if 'G' in code_represents:
                value += 1
            if 'C' in code_represents:
                value += 1

        total_nts += count * 2
        total_gc += count * value

    return (total_gc, total_nts)


def split(sequence, split_by = "123"):
    """Splits a sequence by position, as specified by the 'split_by' parameter. By
    default, the function will split by codon position, and return a dictionary
    containing the keys '1', '2' and '3'.

    The 'split_by' parameter may contain any non-zero number of values, which must
    however be hashable. If a value is specified multiple times, then those positions
    are interleaved (e.g. split_by = "112" returns the first two positions in a codon
    as one sequence, as well as the last positions as one sequence."""
    if not split_by:
        raise TypeError("No partitions to split by specified")

    results = dict((key, []) for key in split_by)
    keys = itertools.chain(itertools.cycle(split_by))
    for (key, nucleotide) in itertools.izip(keys, sequence):
        results[key].append(nucleotide)

    for key in results:
        results[key] = "".join(results[key])

    return results
