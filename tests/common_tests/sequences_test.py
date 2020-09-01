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
import itertools

import pytest

from paleomix.common.sequences import (
    complement,
    reverse_complement,
    encode_genotype,
    split,
)


###############################################################################
###############################################################################
# Tests for 'complement'

_REF_SRC = "ACGTMRWSYKVHDBNX"
_REF_DST = "TGCAKYWSRMBDHVNX"


@pytest.mark.parametrize("src, dst", zip(_REF_SRC, _REF_DST))
def test_complement__single_nt(src, dst):
    assert complement(src) == dst
    assert complement(src.lower()) == dst.lower()


def test_complement__multiple_nts_upper():
    assert complement(_REF_SRC) == _REF_DST


def test_complement__multiple_nts_lower():
    assert complement(_REF_SRC.lower()) == _REF_DST.lower()


def test_complement__multiple_nts_mixed_case():
    assert complement("aGtCn") == "tCaGn"


###############################################################################
###############################################################################
# Tests for 'complement'


def test_reverse_complement():
    assert reverse_complement(_REF_SRC) == _REF_DST[::-1]


###############################################################################
###############################################################################
# Tests for 'encode_genotype'

_IUB_SRC = (
    "A",
    "C",
    "G",
    "T",
    "AC",
    "AG",
    "AT",
    "CG",
    "CT",
    "GT",
    "ACG",
    "ACT",
    "AGT",
    "CGT",
    "ACGT",
)
_IUB_DST = "ACGTMRWSYKVHDB"


@pytest.mark.parametrize("src, dst", zip(_IUB_SRC, _IUB_DST))
def test_genotype__permutations(src, dst):
    for seq in itertools.permutations(src):
        assert encode_genotype("".join(src)) == dst


@pytest.mark.parametrize("value", ("a", "At", "Z", "+"))
def test_genotype__bad_input(value):
    with pytest.raises(ValueError):
        encode_genotype(value)


@pytest.mark.parametrize("sequence", ("CT", "C,T", ",C,T", "C,T,", ",C,T,"))
def test_comma_or_not(sequence):
    assert encode_genotype(sequence) == "Y"


###############################################################################
###############################################################################
# Tests for 'split'


def test_split__empty_sequence():
    assert split("") == {"1": "", "2": "", "3": ""}


def test_split__no_split_by():
    with pytest.raises(ValueError, match="No split_by specified"):
        split("", split_by="")


def test_split__single_group():
    assert split("ACGCAT", "111") == {"1": "ACGCAT"}


def test_split__two_groups():
    assert split("ACGCAT", "112") == {"1": "ACCA", "2": "GT"}


def test_split__three_groups():
    expected = {"1": "AC", "2": "CA", "3": "GT"}
    assert split("ACGCAT", "123") == expected
    assert split("ACGCAT") == expected


def test_split__empty_group():
    expected = {"1": "A", "2": "C", "3": ""}
    assert split("AC") == expected


def test_split__partial_group():
    expected = {"1": "AA", "2": "CA", "3": "G"}
    assert split("ACGAA") == expected
