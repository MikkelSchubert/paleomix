# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import itertools

import pytest

from paleomix.common.sequences import (
    complement,
    encode_genotype,
    reverse_complement,
    split,
)

###############################################################################
###############################################################################
# Tests for 'complement'

_REF_SRC = "ACGTMRWSYKVHDBNX"
_REF_DST = "TGCAKYWSRMBDHVNX"


@pytest.mark.parametrize(("src", "dst"), zip(_REF_SRC, _REF_DST))
def test_complement__single_nt(src: str, dst: str) -> None:
    assert complement(src) == dst
    assert complement(src.lower()) == dst.lower()


def test_complement__multiple_nts_upper() -> None:
    assert complement(_REF_SRC) == _REF_DST


def test_complement__multiple_nts_lower() -> None:
    assert complement(_REF_SRC.lower()) == _REF_DST.lower()


def test_complement__multiple_nts_mixed_case() -> None:
    assert complement("aGtCn") == "tCaGn"


###############################################################################
###############################################################################
# Tests for 'complement'


def test_reverse_complement() -> None:
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


@pytest.mark.parametrize(("src", "dst"), zip(_IUB_SRC, _IUB_DST))
def test_genotype__permutations(src: str, dst: str) -> None:
    for seq in itertools.permutations(src):
        assert encode_genotype("".join(seq)) == dst


@pytest.mark.parametrize("value", ["a", "At", "Z", "-"])
def test_genotype__bad_input(value: str) -> None:
    with pytest.raises(ValueError, match=value):
        encode_genotype(value)


@pytest.mark.parametrize("sequence", ["CT", "C,T", ",C,T", "C,T,", ",C,T,"])
def test_comma_or_not(sequence: str) -> None:
    assert encode_genotype(sequence) == "Y"


###############################################################################
###############################################################################
# Tests for 'split'


def test_split__empty_sequence() -> None:
    assert split("") == {"1": "", "2": "", "3": ""}


def test_split__no_split_by() -> None:
    with pytest.raises(ValueError, match="No split_by specified"):
        split("", split_by="")


def test_split__single_group() -> None:
    assert split("ACGCAT", "111") == {"1": "ACGCAT"}


def test_split__two_groups() -> None:
    assert split("ACGCAT", "112") == {"1": "ACCA", "2": "GT"}


def test_split__three_groups() -> None:
    expected = {"1": "AC", "2": "CA", "3": "GT"}
    assert split("ACGCAT", "123") == expected
    assert split("ACGCAT") == expected


def test_split__empty_group() -> None:
    expected = {"1": "A", "2": "C", "3": ""}
    assert split("AC") == expected


def test_split__partial_group() -> None:
    expected = {"1": "AA", "2": "CA", "3": "G"}
    assert split("ACGAA") == expected
