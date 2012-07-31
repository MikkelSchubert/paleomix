#!/usr/bin/python -3
"""Various functions relating to DNA sequence manipulation."""

from __future__ import print_function

import itertools

 
# Pairs of complementary bases and ambigious basees
_COMPL = [ "AT", "CG", 
           "NN", "RY", 
           "KM", "SS", 
           "WW", "BV", 
           "DH", "XX" ]
_COMPL_TABLE = ["N"] * 256
for (a, b) in _COMPL :
    # Complement both upper/lower-case bases
    for func in (str.upper, str.lower):
        _COMPL_TABLE[ord(func(a))] = func(b)
        _COMPL_TABLE[ord(func(b))] = func(a)
_COMPL_TABLE = "".join(_COMPL_TABLE)



# Table of nt codes used to encode (ambigious) bases
_NT_CODES = [
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
for (abr, nts) in _NT_CODES:
    _NT_CODES_TABLE[frozenset(nts)] = abr
    _NT_CODES_TABLE[frozenset(nts + ",")] = abr

_NT_CODES = dict(_NT_CODES)



def complement(sequence):
    """Returns the complement of a DNA sequence (string)."""
    return sequence.translate(_COMPL_TABLE)


def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence."""
    return complement(sequence)[::-1]


def encode_genotype(nucleotides):
    """Parses a string representing a set of (potentially comma-seperated)
    nucleotides observed at a loci, and returns the corresponding IUPAC code.
    Does not handle lower-case nucleotides, due to lack of clear criteria for
    mixed case input. See e.g. http://www.ebi.ac.uk/2can/tutorials/aa.html"""
    try:
        return _NT_CODES_TABLE[frozenset(nucleotides)]
    except KeyError:
        raise ValueError("Invalid input for 'encode_genotype': '%s'" % (nucleotides, ))



def count_nts(sequence):
    """Given a nucleotide sequence (str), this function returns
    the number of each type of nucleotide representable using 
    IUPAC codes. The sequence must not contain non-IUPAC 
    nucleotides, or other annotation. IUPAC nucleotides are
    handled in a case-insensitive manner."""
    counts = {}
    sequence = sequence.upper()
    for nucleotide in _NT_CODES:
        counts[nucleotide] = sequence.count(nucleotide)

    if len(sequence) != sum(counts.itervalues()):
        raise ValueError("Sequence contains non-(IUPAC-)nucleotides.")
    
    return counts


def count_gc_diploid(sequence):
    """Given a sequence, this function returns a tuple containing the
    the total number of bases that were G/C, as well as the total number
    of bases. The sequence is assumed to represent a diploid genome, with 
    the total number of bases being twice the sequence length, and 
    hence IUPAC codes representing on of or both of G/C are treated as
    reflecting both strands. Thus R counts for 1, while S counts for
    2.

    The sequence must only contain valid IUPAC codes, and no other 
    form of annotation. Both uppercase/lowercase G/Cs are counted."""
    total_nts = total_gc = 0
    counts = count_nts(sequence)
    for (code, count) in counts.iteritems():
        if not count:
            continue

        code_represents = _NT_CODES[code]
        if (len(code_represents) > 2) and (code != 'N'):
            raise ValueError("calculate_gcp assumes diploid genome, nt code for tri-valued SNP observed: " + code)
        
        value = 0
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
