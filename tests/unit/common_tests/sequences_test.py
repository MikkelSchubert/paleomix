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
import unittest
import itertools
import nose.tools
from nose.tools import assert_equal

from pypeline.common.sequences import *



################################################################################
################################################################################
## Tests for 'complement'

_REF_SRC = "ACGTMRWSYKVHDBNX"
_REF_DST = "TGCAKYWSRMBDHVNX"
_REF_PAIRS = zip(_REF_SRC, _REF_DST)


def test_complement__single_nt():
    def test_function(source, destination):
        assert_equal(complement(source), destination)

    for (src, dest) in _REF_PAIRS:
        yield test_function, src, dest
        yield test_function, src.lower(), dest.lower()

def test_complement__multiple_nts_upper():
    assert_equal(complement(_REF_SRC), _REF_DST)

def test_complement__multiple_nts_lower():
    assert_equal(complement(_REF_SRC.lower()), _REF_DST.lower())

def test_complement__multiple_nts_mixed_case():
    assert_equal(complement("aGtCn"), "tCaGn")




################################################################################
################################################################################
## Tests for 'complement'
def test_reverse_complement():
    assert_equal(reverse_complement(_REF_SRC), _REF_DST[::-1])




################################################################################
################################################################################
## Tests for 'genotype'

_IUB_SRC = ("A", "C", "G", "T", 
            "AC", "AG", "AT", "CG", "CT", "GT", 
            "ACG", "ACT", "AGT", "CGT", "ACGT")
_IUB_DST = "ACGTMRWSYKVHDB"
_IUB_PAIRS = zip(_IUB_SRC, _IUB_DST)


def test_genotype__permutations():
    def test_function(src, dst):
        assert_equal(encode_genotype(src), dst)
        
    for (src, dst) in _IUB_PAIRS:
        for seq in itertools.permutations(src):
            yield test_function, "".join(seq), dst


@nose.tools.raises(ValueError)
def test_genotype__bad_input__lowercase():
    encode_genotype("a")

@nose.tools.raises(ValueError)
def test_genotype__bad_input__mixedcase():
    encode_genotype("At")

@nose.tools.raises(ValueError)
def test_genotype__bad_input__unknown_nucleotide():
    encode_genotype("Z")

@nose.tools.raises(ValueError)
def test_genotype__bad_input__non_nucleotide():
    encode_genotype("+")


def test_comma_or_not():
    def test_function(sequence):
        assert_equal(encode_genotype(sequence), "Y")

    for sequence in ("CT", "C,T", ",C,T", "C,T,", ",C,T,"):
        yield test_function, sequence




################################################################################
################################################################################
## Tests for 'split'

def test_split__empty_sequence():
    assert_equal(split(""), {"1" : "", "2" : "", "3" : ""})

@nose.tools.raises(TypeError)
def test_split__no_split_by():
    split("", split_by = "")

def test_split__single_group():
    assert_equal(split("ACGCAT", "111"), {'1' : 'ACGCAT'})

def test_split__two_groups(): 
    assert_equal(split("ACGCAT", "112"), {"1" : "ACCA", "2" : "GT"})

def test_split__three_groups():
    expected = {"1" : "AC", "2" : "CA", "3" : "GT"}
    assert_equal(split("ACGCAT", "123"), expected)
    assert_equal(split("ACGCAT"), expected)

def test_split__empty_group():
    expected = {"1" : "A", "2" : "C", "3" : ""}
    assert_equal(split("AC"), expected)

def test_split__partial_group():
    expected = {"1" : "AA", "2" : "CA", "3" : "G"}
    assert_equal(split("ACGAA"), expected)
