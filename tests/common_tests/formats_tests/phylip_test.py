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
from flexmock import flexmock
from nose.tools import assert_equal

from paleomix.common.formats.phylip import \
     sequential_phy, \
     interleaved_phy

from paleomix.common.formats.msa import \
     MSA

from paleomix.common.formats.fasta import \
     FASTA


_MSA_SHORT_SEQUENCES = \
  MSA([FASTA("seq1", None, "ACGTTGATAACCAGG"),
       FASTA("seq2", None, "TGCAGAGTACGACGT")])
_MSA_MEDIUM_SEQUENCES = \
  MSA([FASTA("seq1", None, "ACGTTGATAACCAGGAGGGATTCGCGATTGGTGGTAACGTAGCC"),
       FASTA("seq2", None, "TGCAGAGTACGACGTCTCCTAGATCCTGGACAATTTAAACCGAA")])
_MSA_LONG_SEQUENCES = \
  MSA([FASTA("seq1", None, "CGGATCTGCTCCTCCACTGGCCACGTTTACTGTCCCCCAACCGTT"
             "CGTCCCGACCTAGTTATACTTCTTAGCAAGGTGTAAAACCAGAGATTGAGGTTATAACG"
             "TTCCTAATCAGTTATTAAATTACCGCGCCCCGACAG"),
       FASTA("seq2", None, "AGTTGAAGAGGCGGAACGTTTGTAAACCGCGCTAACGTAGTTCTA"
             "CAACCAGCCACCCGGTTCGAAGGAACAACTGGTCGCCATAATTAGGCGAAACGATAGTG"
             "CACTAAGGTCAGGTGCGCCCCTGTAAATAATTAGAT")])

_MSA_MEDIUM_NAMES = \
  MSA([FASTA("A_really_long_sequence", None, "ACGTTGATAACCAGG"),
       FASTA("Another_real_long_one!", None, "TGCAGAGTACGACGT")])
_MSA_LONG_NAMES = \
  MSA([FASTA("A_really_long_sequence_name_that_is_in_fact_too_long", None, "ACGTTGATAACCAGG"),
       FASTA("Another_really_long_sequence_name_that_is_too_long", None, "TGCAGAGTACGACGT")])


###############################################################################
###############################################################################
# Tests of 'sequential_phy'

def test_sequential_phy__short_sequences():
    expected = """2 44

seq1
ACGTTGATAA  CCAGGAGGGA  TTCGCGATTG  GTGGTAACGT  AGCC
seq2
TGCAGAGTAC  GACGTCTCCT  AGATCCTGGA  CAATTTAAAC  CGAA"""
    assert_equal(sequential_phy(_MSA_MEDIUM_SEQUENCES), expected)


def test_sequential_phy__multi_line_sequences():
    expected = """2 140

seq1
CGGATCTGCT  CCTCCACTGG  CCACGTTTAC  TGTCCCCCAA  CCGTTCGTCC  CGACCTAGTT
ATACTTCTTA  GCAAGGTGTA  AAACCAGAGA  TTGAGGTTAT  AACGTTCCTA  ATCAGTTATT
AAATTACCGC  GCCCCGACAG
seq2
AGTTGAAGAG  GCGGAACGTT  TGTAAACCGC  GCTAACGTAG  TTCTACAACC  AGCCACCCGG
TTCGAAGGAA  CAACTGGTCG  CCATAATTAG  GCGAAACGAT  AGTGCACTAA  GGTCAGGTGC
GCCCCTGTAA  ATAATTAGAT"""
    assert_equal(sequential_phy(_MSA_LONG_SEQUENCES), expected)


def test_sequential_phy__with_flag():
    expected = """2 15 S

seq1
ACGTTGATAA  CCAGG
seq2
TGCAGAGTAC  GACGT"""
    assert_equal(sequential_phy(_MSA_SHORT_SEQUENCES, add_flag=True), expected)


def test_sequentual_phy__long_names():
    expected = """2 15

A_really_long_sequence_name_th
ACGTTGATAA  CCAGG
Another_really_long_sequence_n
TGCAGAGTAC  GACGT"""
    assert_equal(sequential_phy(_MSA_LONG_NAMES), expected)


def test_sequential_phy__different_lengths():
    _mock = flexmock(MSA).should_receive('validate').at_least.once
    sequential_phy(_MSA_MEDIUM_NAMES)


###############################################################################
###############################################################################
# Tests of 'interleaved_phy'

def test_interleaved_phy__short_sequences():
    expected = """2 44

seq1        ACGTTGATAA  CCAGGAGGGA  TTCGCGATTG  GTGGTAACGT  AGCC
seq2        TGCAGAGTAC  GACGTCTCCT  AGATCCTGGA  CAATTTAAAC  CGAA"""
    assert_equal(interleaved_phy(_MSA_MEDIUM_SEQUENCES), expected)


def test_interleaved_phy__multi_line_sequences():
    expected = """2 140

seq1        CGGATCTGCT  CCTCCACTGG  CCACGTTTAC  TGTCCCCCAA  CCGTTCGTCC
seq2        AGTTGAAGAG  GCGGAACGTT  TGTAAACCGC  GCTAACGTAG  TTCTACAACC

CGACCTAGTT  ATACTTCTTA  GCAAGGTGTA  AAACCAGAGA  TTGAGGTTAT  AACGTTCCTA
AGCCACCCGG  TTCGAAGGAA  CAACTGGTCG  CCATAATTAG  GCGAAACGAT  AGTGCACTAA

ATCAGTTATT  AAATTACCGC  GCCCCGACAG
GGTCAGGTGC  GCCCCTGTAA  ATAATTAGAT"""
    assert_equal(interleaved_phy(_MSA_LONG_SEQUENCES), expected)


def test_interleaved_phy__with_flag():
    expected = """2 15 I

seq1        ACGTTGATAA  CCAGG
seq2        TGCAGAGTAC  GACGT"""
    assert_equal(interleaved_phy(_MSA_SHORT_SEQUENCES, add_flag=True), expected)


def test_interleaved_phy__medium_names():
    expected = """2 15

A_really_long_sequence  ACGTTGATAA  CCAGG
Another_real_long_one!  TGCAGAGTAC  GACGT"""
    assert_equal(interleaved_phy(_MSA_MEDIUM_NAMES), expected)


def test_interleaved_phy__long_names():
    expected = """2 15

A_really_long_sequence_name_th      ACGTTGATAA  CCAGG
Another_really_long_sequence_n      TGCAGAGTAC  GACGT"""
    assert_equal(interleaved_phy(_MSA_LONG_NAMES), expected)


def test_sequentual_phy__different_length_names_1():
    msa = MSA([FASTA("A_short_name", None, "ACGTTGATAACCAGG"),
               FASTA("Another_really_long_sequence_name_that_is_too_long", None, "TGCAGAGTACGACGT")])
    expected = """2 15

A_short_name                        ACGTTGATAA  CCAGG
Another_really_long_sequence_n      TGCAGAGTAC  GACGT"""
    print interleaved_phy(msa), expected
    assert_equal(interleaved_phy(msa), expected)


def test_sequentual_phy__different_length_names_2():
    msa = MSA([FASTA("Burchelli_4", None, "ACGTTGATAACCAGG"),
               FASTA("Donkey",      None, "TGCAGAGTACGACGT")])
    expected = """2 15

Burchelli_4             ACGTTGATAA  CCAGG
Donkey                  TGCAGAGTAC  GACGT"""
    print interleaved_phy(msa), expected
    assert_equal(interleaved_phy(msa), expected)


def test_interleaved_phy__different_lengths():
    _mock = flexmock(MSA).should_receive('validate').at_least.once
    interleaved_phy(_MSA_MEDIUM_NAMES)
