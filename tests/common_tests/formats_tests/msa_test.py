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
import os
import copy
import StringIO

import nose.tools
from nose.tools import \
     assert_equal, \
     assert_raises
from flexmock import \
     flexmock

from paleomix.common.formats.fasta import FASTA
from paleomix.common.formats.msa import \
     MSA, \
     MSAError, \
     FASTAError


def test_dir():
    return os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


def test_file(*args):
    return os.path.join(test_dir(), "data", *args)


###############################################################################
###############################################################################
# Tests for constructor

def test_msa_constructor__calls_validate():
    _mock = flexmock(MSA).should_receive('validate').at_least.once
    MSA([FASTA("NA", None, "ACGT")])


def test_msa_constructor__duplicate_names():
    records = [FASTA("Foo", None, "ACGT"),
               FASTA("Foo", None, "GTCA")]
    assert_raises(MSAError, MSA, records)


def test_msa_constructor__empty_msa():
    assert_raises(MSAError, MSA, [])


###############################################################################
###############################################################################
# Tests for 'seqlen' / len

def test_msa__len__corresponds_to_sequence_number_of_records():
    msa = MSA((FASTA("seq1",    None, "ACGCGTATGCATGCCGA"),
               FASTA("seq2",    None, "TGAACACACAGTAGGAT")))
    assert_equal(len(msa), 2)


def test_msa__seqlen__corresponds_to_sequence_lengths():
    msa = MSA((FASTA("seq1",    None, "ACGCGTATGCATGCCGA"),
               FASTA("seq2",    None, "TGAACACACAGTAGGAT")))
    assert_equal(msa.seqlen(), 17)


###############################################################################
###############################################################################
# Tests for 'exclude'

def test_msa_exclude__remove_one():
    fa_1 = FASTA("A", None, "ACGT")
    fa_2 = FASTA("B", None, "GCTA")
    initial = MSA([fa_1, fa_2])
    expected = MSA([fa_1])
    result = initial.exclude(["B"])
    assert_equal(result, expected)


def test_msa_exclude__missing_keys():
    msa = MSA([FASTA("Foo", None, "ACGT")])
    assert_raises(KeyError, msa.exclude, ["Bar"])


def test_msa_exclude__no_keys():
    msa = MSA([FASTA("Foo", None, "ACGT")])
    assert_raises(ValueError, msa.exclude, [])


###############################################################################
###############################################################################
# Tests for 'select'

def test_msa_select__remove_one():
    fa_1 = FASTA("A", None, "ACGT")
    fa_2 = FASTA("B", None, "GCTA")
    initial = MSA([fa_1, fa_2])
    expected = MSA([fa_1])
    result  = initial.select(["A"])
    assert_equal(result, expected)


def test_msa_select__missing_keys():
    msa = MSA([FASTA("Foo", None, "ACGT")])
    assert_raises(KeyError, msa.select, ["Bar"])


def test_msa_select__no_keys():
    msa = MSA([FASTA("Foo", None, "ACGT")])
    assert_raises(ValueError, msa.select, [])


###############################################################################
###############################################################################
# Tests for 'reduce'

def test_msa_reduce__no_empty_columns__no_columns_are_removed():
    fa_1 = FASTA("Name_A", "Meta_A", "ACnT")
    fa_2 = FASTA("Name_B", "Meta_B", "C-TN")
    initial = MSA([fa_1, fa_2])
    expected = MSA([fa_1, fa_2])
    assert_equal(initial.reduce(), expected)


def test_msa_reduce__one_empty_column__column_are_removed():
    fa_1 = FASTA("Name_A", "Meta_A", "AnT")
    fa_2 = FASTA("Name_B", "Meta_B", "C-N")
    initial = MSA([fa_1, fa_2])
    fa_reduced_1 = FASTA("Name_A", "Meta_A", "AT")
    fa_reduced_2 = FASTA("Name_B", "Meta_B", "CN")
    expected = MSA([fa_reduced_1, fa_reduced_2])
    assert_equal(initial.reduce(), expected)


def test_msa_reduce__multiple_empty_column__all_empty_column_are_removed():
    fa_1 = FASTA("Name_A", "Meta_A", "-AnTN")
    fa_2 = FASTA("Name_B", "Meta_B", "NC-NN")
    initial = MSA([fa_1, fa_2])
    fa_reduced_1 = FASTA("Name_A", "Meta_A", "AT")
    fa_reduced_2 = FASTA("Name_B", "Meta_B", "CN")
    expected = MSA([fa_reduced_1, fa_reduced_2])
    assert_equal(initial.reduce(), expected)


def test_msa_reduce__only_empty_column__none_is_returned():
    fa_1 = FASTA("Name_A", "Meta_A", "---Nn")
    fa_2 = FASTA("Name_B", "Meta_B", "Nn--N")
    initial = MSA([fa_1, fa_2])
    assert_equal(initial.reduce(), None)


###############################################################################
###############################################################################
# Tests for 'filter_singletons'

_FILTER_MSA_1 = MSA((FASTA("Seq1", "Meta1", "ACGNTYCSTG"),
                     FASTA("Seq2", "Meta2", "ACTA-WCCTG"),
                     FASTA("Seq3", "Meta3", "NCGGTYCGTC")))


def test_msa_filter_singletons__filter_by_second():
    expected = MSA((FASTA("Seq1", "Meta1", "ACnNntCcTG"),
                    FASTA("Seq2", "Meta2", "ACTA-WCCTG"),
                    FASTA("Seq3", "Meta3", "NCGGTYCGTC")))
    result = _FILTER_MSA_1.filter_singletons("Seq1", ["Seq2"])
    assert_equal(result, expected)


def test_msa_filter_singletons__filter_by_third():
    expected = MSA((FASTA("Seq1", "Meta1", "nCGNTYCgTn"),
                    FASTA("Seq2", "Meta2", "ACTA-WCCTG"),
                    FASTA("Seq3", "Meta3", "NCGGTYCGTC")))
    result = _FILTER_MSA_1.filter_singletons("Seq1", ["Seq3"])
    assert_equal(result, expected)


def test_msa_filter_singletons__filter_by_both():
    result = _FILTER_MSA_1.filter_singletons("Seq1", ["Seq2", "Seq3"])
    assert_equal(result, _FILTER_MSA_1)


def test_msa_filter_singletons__filter_by_itself():
    assert_raises(MSAError, _FILTER_MSA_1.filter_singletons, "Seq1", ["Seq1", "Seq2"])


def test_msa_filter_singletons__filter_by_nothing():
    assert_raises(ValueError, _FILTER_MSA_1.filter_singletons, "Seq1", [])


###############################################################################
###############################################################################
# Tests for 'MSA.join'

_JOIN_MSA_1 = MSA((FASTA("nc",    None, "ACG"),
                   FASTA("nm",    None, "TGA"),
                   FASTA("miRNA", None, "UCA")))
_JOIN_MSA_2 = MSA((FASTA("nc",    None, "TGA"),
                   FASTA("nm",    None, "CTT"),
                   FASTA("miRNA", None, "GAC")))
_JOIN_MSA_3 = MSA((FASTA("nc",    None, "AAG"),
                   FASTA("nm",    None, "GAG"),
                   FASTA("miRNA", None, "CAU")))


def test_msa_join__single_msa():
    result = MSA.join(_JOIN_MSA_1)
    assert_equal(result, _JOIN_MSA_1)


def test_msa_join__two_msa():
    expected = MSA((FASTA("nc",    None, "ACGTGA"),
                    FASTA("nm",    None, "TGACTT"),
                    FASTA("miRNA", None, "UCAGAC")))
    result = MSA.join(_JOIN_MSA_1, _JOIN_MSA_2)
    assert_equal(result, expected)


def test_msa_join__three_msa():
    expected = MSA((FASTA("nc",    None, "ACGTGAAAG"),
                    FASTA("nm",    None, "TGACTTGAG"),
                    FASTA("miRNA", None, "UCAGACCAU")))
    result = MSA.join(_JOIN_MSA_1, _JOIN_MSA_2, _JOIN_MSA_3)
    assert_equal(result, expected)


@nose.tools.raises(TypeError)
def test_msa_join__missing_arguments():
    MSA.join()


###############################################################################
###############################################################################
# Tests for 'MSA.from_lines'

def test_msa_from_lines__single_entry():
    lines = [">seq1", "ACG"]
    result = MSA([FASTA("seq1", None, "ACG")])
    assert_equal(MSA.from_lines(lines), result)


def test_msa_from_lines__single_entry_with_meta():
    lines = [">seq1 Meta info", "ACG"]
    expected = MSA([FASTA("seq1", "Meta info", "ACG")])
    result = MSA.from_lines(lines)
    assert_equal(result, expected)


def test_msa_from_lines__two_entries():
    lines = [">seq1", "ACG", ">seq2", "TGA"]
    expected = MSA([FASTA("seq1", None, "ACG"),
                    FASTA("seq2", None, "TGA")])
    result = MSA.from_lines(lines)
    assert_equal(result, expected)


def test_msa_from_lines__two_entries_with_meta():
    lines = [">seq1", "ACG", ">seq2 Second meta", "TGA"]
    expected = MSA([FASTA("seq1", None, "ACG"),
                    FASTA("seq2", "Second meta", "TGA")])
    result = MSA.from_lines(lines)
    assert_equal(result, expected)


@nose.tools.raises(MSAError)
def test_msa_from_lines__duplicate_names():
    MSA.from_lines([">seq1", "ACG", ">seq1", "TGA"])


@nose.tools.raises(MSAError)
def test_msa_from_lines__mismatched_lengths():
    MSA.from_lines([">seq1", "ACG", ">seq2", "TGAN"])


@nose.tools.raises(FASTAError)
def test_msa_from_lines__empty_name():
    MSA.from_lines([">", "ACG", ">seq1", "TGAN"])


###############################################################################
###############################################################################
# Tests for 'MSA.from_file'

def test_msa_from_file__uncompressed():
    expected = MSA([FASTA("This_is_FASTA!", None, "ACGTN"),
                    FASTA("This_is_ALSO_FASTA!", None, "CGTNA")])
    results = MSA.from_file(test_file("fasta_file.fasta"))
    assert_equal(results, expected)


def test_msa_from_file__compressed_gz():
    expected = MSA([FASTA("This_is_GZipped_FASTA!", None, "ACGTN"),
                    FASTA("This_is_ALSO_GZipped_FASTA!", None, "CGTNA")])
    results = MSA.from_file(test_file("fasta_file.fasta.gz"))
    assert_equal(results, expected)


def test_msa_from_file__compressed_bz2():
    expected = MSA([FASTA("This_is_BZ_FASTA!", None, "CGTNA"),
                    FASTA("This_is_ALSO_BZ_FASTA!", None,  "ACGTN")])
    results = MSA.from_file(test_file("fasta_file.fasta.bz2"))
    assert_equal(results, expected)


###############################################################################
###############################################################################
# Tests for 'MSA.split'

def test_msa_split_msa__single_group():
    msa = MSA([FASTA("seq1", None, "ACGCAT"),
               FASTA("seq2", None, "GAGTGA")])
    expected = {'1': copy.copy(msa)}
    assert_equal(msa.split("111"), expected)


def test_msa_split_msa__two_groups():
    msa = MSA([FASTA("seq1", None, "ACGCAT"),
               FASTA("seq2", None, "GAGTGA")])
    expected = {"1": MSA([FASTA("seq1", None, "ACCA"),
                          FASTA("seq2", None, "GATG")]),
                "2": MSA([FASTA("seq1", None, "GT"),
                          FASTA("seq2", None, "GA")])}
    assert_equal(msa.split("112"), expected)


def test_msa_split__three_groups():
    msa = MSA([FASTA("seq1", None, "ACGCAT"),
               FASTA("seq2", None, "GAGTGA")])
    expected = {"1": MSA([FASTA("seq1", None, "AC"),
                          FASTA("seq2", None, "GT")]),
                "2": MSA([FASTA("seq1", None, "CA"),
                          FASTA("seq2", None, "AG")]),
                "3": MSA([FASTA("seq1", None, "GT"),
                          FASTA("seq2", None, "GA")])}
    assert_equal(msa.split("123"), expected)


def test_msa_split__empty_group():
    msa = MSA([FASTA("seq1", None, "AC"),
               FASTA("seq2", None, "GA")])
    expected = {"1": MSA([FASTA("seq1", None, "A"),
                          FASTA("seq2", None, "G")]),
                "2": MSA([FASTA("seq1", None, "C"),
                          FASTA("seq2", None, "A")]),
                "3": MSA([FASTA("seq1", None, ""),
                          FASTA("seq2", None, "")])}
    assert_equal(msa.split("123"), expected)


def test_msa_split__partial_group():
    msa = MSA([FASTA("seq1", None, "ACGCA"),
               FASTA("seq2", None, "GAGTG")])
    expected = {"1": MSA([FASTA("seq1", None, "AC"),
                          FASTA("seq2", None, "GT")]),
                "2": MSA([FASTA("seq1", None, "CA"),
                          FASTA("seq2", None, "AG")]),
                "3": MSA([FASTA("seq1", None, "G"),
                          FASTA("seq2", None, "G")])}
    assert_equal(msa.split("123"), expected)


@nose.tools.raises(TypeError)
def test_msa_split_msa__no_split_by():
    msa = MSA([FASTA("seq1", None, "ACG"),
               FASTA("seq2", None, "GAT")])
    msa.split(split_by="")


###############################################################################
###############################################################################
# Tests for 'MSA.to_file'

def test_msa_to_file__complete_line_test():
    msa = MSA([FASTA("barfoo", None, "ACGATA" * 10 + "CGATAG" * 5),
               FASTA("foobar", None, "CGAATG" * 10 + "TGTCAT" * 5)])
    expected = ">barfoo\n%s\n%s\n" % ("ACGATA" * 10, "CGATAG" * 5)
    expected += ">foobar\n%s\n%s\n" % ("CGAATG" * 10, "TGTCAT" * 5)
    stringf = StringIO.StringIO()
    MSA.to_file(msa, stringf)
    assert_equal(stringf.getvalue(), expected)


###############################################################################
###############################################################################
# Tests for 'MSA.validate'

def test_msa_validate__missing_names_first():
    msa_1 = MSA(list(_JOIN_MSA_1)[:-1])
    msa_2 = copy.copy(_JOIN_MSA_2)
    assert_raises(MSAError, MSA.validate, msa_1, msa_2)


def test_msa_validate__missing_names_second():
    msa_1 = copy.copy(_JOIN_MSA_1)
    msa_2 = MSA(list(_JOIN_MSA_2)[:-1])
    assert_raises(MSAError, MSA.validate, msa_1, msa_2)


###############################################################################
###############################################################################
# Tests for 'MSA.names'

def test_msa_names():
    assert_equal(_JOIN_MSA_1.names(), set(("nc", "nm", "miRNA")))


###############################################################################
###############################################################################
# Tests for str/repr

def test_msa_repr():
    msa = MSA((FASTA("nc",    None, "ACGTA"),
               FASTA("nm",    "META", "TGAGT"),
               FASTA("miRNA", None, "UCAGA")))
    expected = "MSA(FASTA('miRNA', None, 'UCAGA'), FASTA('nc', None, 'ACGTA'), FASTA('nm', 'META', 'TGAGT'))"
    assert_equal(str(msa), expected)


def test_msa_repr__same_as_str():
    assert_equal(str(_JOIN_MSA_1), repr(_JOIN_MSA_1))
