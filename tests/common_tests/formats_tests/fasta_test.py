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
import StringIO

import nose.tools
from nose.tools import \
    assert_is, \
    assert_equal, \
    assert_raises, \
    assert_not_equal, \
    assert_less, \
    assert_less_equal, \
    assert_greater, \
    assert_greater_equal


from paleomix.common.testing import assert_list_equal
from paleomix.common.formats.fasta import \
     FASTA, \
     FASTAError


###############################################################################
###############################################################################

_SEQ_FRAG = "AAGTCC"  # len() = 6


def test_dir():
    return os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


def test_file(*args):
    return os.path.join(test_dir(), "data", *args)


###############################################################################
###############################################################################
# Tests for FASTA constructor

def _simple_fasta_record():
    return FASTA("Dummy", "Meta-inf", "ACGT")


def test_fasta__constructor__name():
    record = _simple_fasta_record()
    assert_equal(record.name, "Dummy")


def test_fasta__constructor__meta():
    record = _simple_fasta_record()
    assert_equal(record.meta, "Meta-inf")


def test_fasta__constructor__sequence():
    record = _simple_fasta_record()
    assert_equal(record.sequence, "ACGT")


def test_fasta__constructor__name_must_be_non_empty():
    assert_raises(FASTAError, FASTA, "", None, "ACGT")


def test_fasta__constructor__name_must_be_string_type():
    assert_raises(FASTAError, FASTA, 1, None, "ACGT")


def test_fasta__constructor__name_must_be_string_type_or_none():
    assert_raises(FASTAError, FASTA, "Seq1", 1, "ACGT")


def test_fasta__constructor__sequence_must_be_string_type():
    assert_raises(FASTAError, FASTA, "Seq1", None, 1)


###############################################################################
###############################################################################
# Tests for __repr__

def test_fasta__repr__partial_line_test():
    expected = ">foobar\n%s\n" % (_SEQ_FRAG, )
    result = repr(FASTA("foobar", None, _SEQ_FRAG))
    assert_equal(result, expected)


def test_fasta__repr__complete_line_test():
    expected = ">barfoo\n%s\n" % (_SEQ_FRAG * 10, )
    result = repr(FASTA("barfoo", None, _SEQ_FRAG * 10))
    assert_equal(result, expected)


def test_fasta__repr__multiple_lines():
    expected = ">foobar\n%s\n%s\n" \
        % (_SEQ_FRAG * 10, _SEQ_FRAG * 5)
    result = repr(FASTA("foobar", None, _SEQ_FRAG * 15))
    assert_equal(result, expected)


def test_fasta__repr__partial_line_test_with_meta_information():
    expected = ">foobar my Meta-Info\n%s\n" % (_SEQ_FRAG, )
    result = repr(FASTA("foobar", "my Meta-Info", _SEQ_FRAG))
    assert_equal(result, expected)


###############################################################################
###############################################################################
# Tests for print_fasta


def test_fasta__write__partial_line():
    expected = ">foobar\n%s\n" % (_SEQ_FRAG, )
    stringf = StringIO.StringIO()
    FASTA("foobar", None, _SEQ_FRAG).write(stringf)
    assert_equal(stringf.getvalue(), expected)


def test_fasta__write__complete_line_test():
    expected = ">barfoo\n%s\n" % (_SEQ_FRAG * 10, )
    stringf = StringIO.StringIO()
    FASTA("barfoo", None, _SEQ_FRAG * 10).write(stringf)
    assert_equal(stringf.getvalue(), expected)


def test_fasta__write__multiple_lines():
    expected = ">foobar\n%s\n%s\n" \
        % (_SEQ_FRAG * 10, _SEQ_FRAG * 5)
    stringf = StringIO.StringIO()
    FASTA("foobar", None, _SEQ_FRAG * 15).write(stringf)
    assert_equal(stringf.getvalue(), expected)


###############################################################################
###############################################################################
# Tests for FASTA.from_lines

def test_fasta__from_lines__no_records():
    assert_list_equal(FASTA.from_lines([]), [])


def test_fasta__from_lines_single_record():
    lines = [">single\n", "TGTTCTCCACCGTGCACAAC\n", "CCTTCATCCA\n"]
    expected = [FASTA("single", None, "TGTTCTCCACCGTGCACAACCCTTCATCCA")]
    assert_list_equal(FASTA.from_lines(lines), expected)


def test_fasta__from_lines__multiple_records():
    lines = [">first\n",  "TGTTCTCCACCGTGCACAAC\n", "CCTTCATCCA\n",
             ">Second XT:1:0\n", "GAGAGCTCAGCTAAC\n",
             ">Third\n",  "CGCTGACCAAAAACGGACAG\n", "GGCATTCGGC\n"]
    expected = [FASTA("first", None, "TGTTCTCCACCGTGCACAACCCTTCATCCA"),
                FASTA("Second", "XT:1:0", "GAGAGCTCAGCTAAC"),
                FASTA("Third", None, "CGCTGACCAAAAACGGACAGGGCATTCGGC")]
    assert_list_equal(FASTA.from_lines(lines), expected)


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__empty_record_name_only__nothing_else():
    list(FASTA.from_lines([">fasta1\n"]))


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__empty_record_name_only__first():
    list(FASTA.from_lines([">fasta1\n", ">fasta2\n", "AGTC\n"]))


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__empty_record__middle():
    lines = [">fasta0\n", "ACGT\n", ">fasta1\n", ">fasta2\n", "AGTC\n"]
    list(FASTA.from_lines(lines))


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__empty_record_last():
    lines = [">fasta1\n", "ACGT\n", ">fasta2\n"]
    list(FASTA.from_lines(lines))


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__missing_name__alone():
    lines = ["ACGT\n"]
    list(FASTA.from_lines(lines))


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__missing_name__with_others():
    lines = ["ACGT\n", ">Foo\n", "ACGGTA\n"]
    list(FASTA.from_lines(lines))


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__empty_name__alone():
    lines = [">\n", "ACGT\n"]
    list(FASTA.from_lines(lines))


@nose.tools.raises(FASTAError)
def test_fasta__from_lines__empty_name__with_others():
    lines = [">\n", "ACGT\n", ">Foo\n", "ACGGTA\n"]
    list(FASTA.from_lines(lines))


###############################################################################
###############################################################################
# Tests for 'FASTA.from_file'

def test_fasta__from_file__uncompressed():
    expected = [FASTA("This_is_FASTA!", None, "ACGTN"),
                FASTA("This_is_ALSO_FASTA!", None, "CGTNA")]
    results = list(FASTA.from_file(test_file("fasta_file.fasta")))
    assert_equal(results, expected)


def test_fasta__from_file__compressed_gz():
    expected = [FASTA("This_is_GZipped_FASTA!", None, "ACGTN"),
                FASTA("This_is_ALSO_GZipped_FASTA!", None, "CGTNA")]
    results = list(FASTA.from_file(test_file("fasta_file.fasta.gz")))
    assert_equal(results, expected)


def test_fasta__from_file__compressed_bz2():
    expected = [FASTA("This_is_BZ_FASTA!", None, "CGTNA"),
                FASTA("This_is_ALSO_BZ_FASTA!", None, "ACGTN")]
    results = list(FASTA.from_file(test_file("fasta_file.fasta.bz2")))
    assert_equal(results, expected)


###############################################################################
###############################################################################

def test_fasta__equality():
    assert_equal(FASTA("A", "B", "C"), FASTA("A", "B", "C"))


def test_fasta__inequality():
    assert_not_equal(FASTA("A", "B", "C"), FASTA("A", "B", "D"))
    assert_not_equal(FASTA("A", "B", "C"), FASTA("A", None, "C"))
    assert_not_equal(FASTA("A", "B", "C"), FASTA("D", "B", "C"))


def test_fasta__sorting_less_equal():
    assert not FASTA("A", "B", "C") < FASTA("A", "B", "C")
    assert_less(FASTA("A", "B", "C"), FASTA("B", "B", "C"))
    assert_less(FASTA("A", "B", "C"), FASTA("A", "C", "C"))
    assert_less(FASTA("A", "B", "C"), FASTA("A", "B", "D"))
    assert_less_equal(FASTA("A", "B", "C"), FASTA("A", "B", "C"))
    assert_less_equal(FASTA("A", "B", "C"), FASTA("B", "B", "C"))
    assert_less_equal(FASTA("A", "B", "C"), FASTA("A", "C", "C"))
    assert_less_equal(FASTA("A", "B", "C"), FASTA("A", "B", "D"))


def test_fasta__sorting_greater_equal():
    assert not FASTA("A", "B", "C") > FASTA("A", "B", "C")
    assert_greater(FASTA("B", "B", "C"), FASTA("A", "B", "C"))
    assert_greater(FASTA("A", "C", "C"), FASTA("A", "B", "C"))
    assert_greater(FASTA("A", "B", "D"), FASTA("A", "B", "C"))
    assert_greater_equal(FASTA("A", "B", "C"), FASTA("A", "B", "C"))
    assert_greater_equal(FASTA("B", "B", "C"), FASTA("A", "B", "C"))
    assert_greater_equal(FASTA("A", "C", "C"), FASTA("A", "B", "C"))
    assert_greater_equal(FASTA("A", "B", "D"), FASTA("A", "B", "C"))


def test_fasta__hash():
    assert_equal(hash(FASTA("A", "B", "C")), hash(FASTA("A", "B", "C")))
    assert_not_equal(hash(FASTA("A", "B", "C")), hash(FASTA("B", "B", "C")))
    assert_not_equal(hash(FASTA("A", "B", "C")), hash(FASTA("A", "C", "C")))
    assert_not_equal(hash(FASTA("A", "B", "C")), hash(FASTA("A", "B", "D")))


def test_fasta__unimplemented_comparison():
    assert_is(NotImplemented, FASTA("A", None, "C").__eq__(10))
    assert_is(NotImplemented, FASTA("A", None, "C").__lt__(10))
    assert_is(NotImplemented, FASTA("A", None, "C").__le__(10))
    assert_is(NotImplemented, FASTA("A", None, "C").__ge__(10))
    assert_is(NotImplemented, FASTA("A", None, "C").__gt__(10))

