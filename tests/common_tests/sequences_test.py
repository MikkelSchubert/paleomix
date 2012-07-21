import unittest
import itertools
import nose.tools

from common.sequences import *



################################################################################
################################################################################
## Tests for 'complement'

_REF_SRC = "ACGTMRWSYKVHDBNX"
_REF_DST = "TGCAKYWSRMBDHVNX"
_REF_PAIRS = zip(_REF_SRC, _REF_DST)


def test_complement__single_nt():
    def test_function(source, destination):
        assert complement(source) == destination

    for (src, dest) in _REF_PAIRS:
        yield test_function, src, dest
        yield test_function, src.lower(), dest.lower()

def test_complement__multiple_nts_upper():
    assert complement(_REF_SRC) == _REF_DST

def test_complement__multiple_nts_lower():
    assert complement(_REF_SRC.lower()) == _REF_DST.lower()

def test_complement__multiple_nts_mixed_case():
    assert complement("aGtCn") == "tCaGn"


################################################################################
################################################################################
## Tests for 'complement'
def test_reverse_complement():
    assert reverse_complement(_REF_SRC) == _REF_DST[::-1]


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
        assert encode_genotype(src) == dst
        
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
        assert encode_genotype(sequence) == "Y"

    for sequence in ("CT", "C,T", ",C,T", "C,T,", ",C,T,"):
        yield test_function, sequence
