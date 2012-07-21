#!/usr/bin/python -3
"""Various functions relating to DNA sequence manipulation."""

from __future__ import print_function

import sys
import gzip
import string

from pypeline.common.utilities import grouper, split_before
from pypeline.common.formats.fasta import *

 
# Pairs of complementary bases and ambigious basees
_compl = [ "AT", "CG", 
           "NN", "RY", 
           "KM", "SS", 
           "WW", "BV", 
           "DH", "XX" ]
_compl_table = ["N"] * 256
for (a, b) in _compl:
    # Complement both upper/lower-case bases
    for func in (str.upper, str.lower):
        _compl_table[ord(func(a))] = func(b)
        _compl_table[ord(func(b))] = func(a)
_compl_table = "".join(_compl_table)



# Table of nt codes used to encode (ambigious) bases
_nt_codes = [
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

_nt_codes_table = {}
for (abr, nts) in _nt_codes:
    _nt_codes_table[frozenset(nts)] = abr
    _nt_codes_table[frozenset(nts + ",")] = abr

_nt_codes = dict(_nt_codes)



def complement(s):
    """Returns the complement of a DNA sequence (string)."""
    return s.translate(_compl_table)


def reverse_complement(s):
    """Returns the reverse complement of a DNA sequence."""
    return complement(s)[::-1]


def encode_genotype(nts):
    """Parses a string representing a set of (potentially comma-seperated)
    nucleotides observed at a loci, and returns the corresponding IUPAC code.
    Does not handle lower-case nucleotides, due to lack of clear criteria for
    mixed case input. See e.g. http://www.ebi.ac.uk/2can/tutorials/aa.html"""
    try:
        return _nt_codes_table[frozenset(nts)]
    except KeyError:
        raise ValueError("Invalid input for 'encode_genotype': '%s'" % (nts, ))



def count_nts(sequence):
    """Given a nucleotide sequence (str), this function returns
    the number of each type of nucleotide representable using 
    IUPAC codes. The sequence must not contain non-IUPAC 
    nucleotides, or other annotation. IUPAC nucleotides are
    handled in a case-insensitive manner."""
    counts = {}
    sequence = sequence.upper()
    for nt in _nt_codes:
        counts[nt] = sequence.count(nt)

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

        code_represents = _nt_codes[code]
        if (len(code_represents) > 2) and (code != 'N'):
            raise ValueError("calculate_gcp assumes diploid genome, nt code for tri-valued SNP observed: " + code)
        
        value = 0
        if 'G' in code_represents: value += 1
        if 'C' in code_represents: value += 1

        total_nts += count * 2
        total_gc += count * value

    return (total_gc, total_nts)
