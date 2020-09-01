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
import bz2
import gzip
import io

import pytest

from paleomix.common.fileutils import fspath
from paleomix.common.formats.fasta import FASTA, FASTAError


_SEQ_FRAG = "AAGTCC"  # len() = 6


###############################################################################
###############################################################################
# Tests for FASTA constructor


def _simple_fasta_record():
    return FASTA("Dummy", "Meta-inf", "ACGT")


def test_fasta__constructor__name():
    record = _simple_fasta_record()
    assert record.name == "Dummy"


def test_fasta__constructor__meta():
    record = _simple_fasta_record()
    assert record.meta == "Meta-inf"


def test_fasta__constructor__sequence():
    record = _simple_fasta_record()
    assert record.sequence == "ACGT"


def test_fasta__constructor__name_must_be_non_empty():
    with pytest.raises(FASTAError):
        FASTA("", None, "ACGT")


def test_fasta__constructor__name_must_be_string_type():
    with pytest.raises(FASTAError):
        FASTA(1, None, "ACGT")


def test_fasta__constructor__name_must_be_string_type_or_none():
    with pytest.raises(FASTAError):
        FASTA("Seq1", 1, "ACGT")


def test_fasta__constructor__sequence_must_be_string_type():
    with pytest.raises(FASTAError):
        FASTA("Seq1", None, 1)


###############################################################################
###############################################################################
# Tests for __repr__


def test_fasta__repr__partial_line_test():
    expected = "FASTA('foobar', '', %r)" % (_SEQ_FRAG,)
    result = repr(FASTA("foobar", None, _SEQ_FRAG))
    assert result == expected


def test_fasta__repr__complete_line_test():
    expected = "FASTA('barfoo', '', %r)" % (_SEQ_FRAG * 10,)
    result = repr(FASTA("barfoo", None, _SEQ_FRAG * 10))
    assert result == expected


def test_fasta__repr__multiple_lines():
    expected = "FASTA('foobar', '', %r)" % (_SEQ_FRAG * 15,)
    result = repr(FASTA("foobar", None, _SEQ_FRAG * 15))
    assert result == expected


def test_fasta__repr__partial_line_test_with_meta_information():
    expected = "FASTA('foobar', 'my Meta-Info', %r)" % (_SEQ_FRAG,)
    result = repr(FASTA("foobar", "my Meta-Info", _SEQ_FRAG))
    assert result == expected


###############################################################################
###############################################################################
# Tests for write


def test_fasta__write__partial_line():
    expected = ">foobar\n%s\n" % (_SEQ_FRAG,)
    stringf = io.StringIO()
    FASTA("foobar", None, _SEQ_FRAG).write(stringf)
    assert stringf.getvalue() == expected


def test_fasta__write__with_metadata():
    expected = ">foobar my Meta data\n%s\n" % (_SEQ_FRAG,)
    stringf = io.StringIO()
    FASTA("foobar", "my Meta data", _SEQ_FRAG).write(stringf)
    assert stringf.getvalue() == expected


def test_fasta__write__complete_line_test():
    expected = ">barfoo\n%s\n" % (_SEQ_FRAG * 10,)
    stringf = io.StringIO()
    FASTA("barfoo", None, _SEQ_FRAG * 10).write(stringf)
    assert stringf.getvalue() == expected


def test_fasta__write__multiple_lines():
    expected = ">foobar\n%s\n%s\n" % (_SEQ_FRAG * 10, _SEQ_FRAG * 5)
    stringf = io.StringIO()
    FASTA("foobar", None, _SEQ_FRAG * 15).write(stringf)
    assert stringf.getvalue() == expected


###############################################################################
###############################################################################
# Tests for FASTA.from_lines


def test_fasta__from_lines__no_records():
    assert list(FASTA.from_lines([])) == list([])


def test_fasta__from_lines_single_record():
    lines = [">single\n", "TGTTCTCCACCGTGCACAAC\n", "CCTTCATCCA\n"]
    expected = [FASTA("single", None, "TGTTCTCCACCGTGCACAACCCTTCATCCA")]
    assert list(FASTA.from_lines(lines)) == list(expected)


def test_fasta__from_lines__multiple_records():
    lines = [
        ">first\n",
        "TGTTCTCCACCGTGCACAAC\n",
        "CCTTCATCCA\n",
        ">Second XT:1:0\n",
        "GAGAGCTCAGCTAAC\n",
        ">Third\n",
        "CGCTGACCAAAAACGGACAG\n",
        "GGCATTCGGC\n",
    ]
    expected = [
        FASTA("first", None, "TGTTCTCCACCGTGCACAACCCTTCATCCA"),
        FASTA("Second", "XT:1:0", "GAGAGCTCAGCTAAC"),
        FASTA("Third", None, "CGCTGACCAAAAACGGACAGGGCATTCGGC"),
    ]
    assert list(FASTA.from_lines(lines)) == list(expected)


def test_fasta__from_lines__empty_record_name_only__nothing_else():
    with pytest.raises(FASTAError):
        list(FASTA.from_lines([">fasta1\n"]))


def test_fasta__from_lines__empty_record_name_only__first():
    with pytest.raises(FASTAError):
        list(FASTA.from_lines([">fasta1\n", ">fasta2\n", "AGTC\n"]))


def test_fasta__from_lines__empty_record__middle():
    lines = [">fasta0\n", "ACGT\n", ">fasta1\n", ">fasta2\n", "AGTC\n"]
    with pytest.raises(FASTAError):
        list(FASTA.from_lines(lines))


def test_fasta__from_lines__empty_record_last():
    lines = [">fasta1\n", "ACGT\n", ">fasta2\n"]
    with pytest.raises(FASTAError):
        list(FASTA.from_lines(lines))


def test_fasta__from_lines__missing_name__alone():
    with pytest.raises(FASTAError):
        list(FASTA.from_lines(["ACGT\n"]))


def test_fasta__from_lines__missing_name__with_others():
    with pytest.raises(FASTAError):
        list(FASTA.from_lines(["ACGT\n", ">Foo\n", "ACGGTA\n"]))


def test_fasta__from_lines__empty_name__alone():
    with pytest.raises(FASTAError):
        list(FASTA.from_lines([">\n", "ACGT\n"]))


def test_fasta__from_lines__empty_name__with_others():
    with pytest.raises(FASTAError):
        list(FASTA.from_lines([">\n", "ACGT\n", ">Foo\n", "ACGGTA\n"]))


###############################################################################
###############################################################################
# Tests for 'FASTA.from_file'


@pytest.mark.parametrize("func", (open, gzip.open, bz2.open))
def test_fasta__from_file(func, tmp_path):
    expected = [
        FASTA("This_is_FASTA!", None, "ACGTN"),
        FASTA("This_is_ALSO_FASTA!", None, "CGTNA"),
    ]

    with func(fspath(tmp_path / "file"), "wt") as handle:
        for item in expected:
            item.write(handle)

    assert list(FASTA.from_file(tmp_path / "file")) == expected


###############################################################################
###############################################################################


def test_fasta__equality():
    assert FASTA("A", "B", "C") == FASTA("A", "B", "C")


def test_fasta__inequality():
    assert FASTA("A", "B", "C") != FASTA("A", "B", "D")
    assert FASTA("A", "B", "C") != FASTA("A", None, "C")
    assert FASTA("A", "B", "C") != FASTA("D", "B", "C")


def test_fasta__sorting_less_equal():
    assert not FASTA("A", "B", "C") < FASTA("A", "B", "C")
    assert FASTA("A", "B", "C") < FASTA("B", "B", "C")
    assert FASTA("A", "B", "C") < FASTA("A", "C", "C")
    assert FASTA("A", "B", "C") < FASTA("A", "B", "D")
    assert FASTA("A", "B", "C") <= FASTA("A", "B", "C")
    assert FASTA("A", "B", "C") <= FASTA("B", "B", "C")
    assert FASTA("A", "B", "C") <= FASTA("A", "C", "C")
    assert FASTA("A", "B", "C") <= FASTA("A", "B", "D")


def test_fasta__sorting_greater_equal():
    assert not FASTA("A", "B", "C") > FASTA("A", "B", "C")
    assert FASTA("B", "B", "C") > FASTA("A", "B", "C")
    assert FASTA("A", "C", "C") > FASTA("A", "B", "C")
    assert FASTA("A", "B", "D") > FASTA("A", "B", "C")
    assert FASTA("A", "B", "C") >= FASTA("A", "B", "C")
    assert FASTA("B", "B", "C") >= FASTA("A", "B", "C")
    assert FASTA("A", "C", "C") >= FASTA("A", "B", "C")
    assert FASTA("A", "B", "D") >= FASTA("A", "B", "C")


def test_fasta__hash():
    assert hash(FASTA("A", "B", "C")) == hash(FASTA("A", "B", "C"))
    assert hash(FASTA("A", "B", "C")) != hash(FASTA("B", "B", "C"))
    assert hash(FASTA("A", "B", "C")) != hash(FASTA("A", "C", "C"))
    assert hash(FASTA("A", "B", "C")) != hash(FASTA("A", "B", "D"))


def test_fasta__unimplemented_comparison():
    assert NotImplemented is FASTA("A", None, "C").__eq__(10)
    assert NotImplemented is FASTA("A", None, "C").__lt__(10)
    assert NotImplemented is FASTA("A", None, "C").__le__(10)
    assert NotImplemented is FASTA("A", None, "C").__ge__(10)
    assert NotImplemented is FASTA("A", None, "C").__gt__(10)


###############################################################################
###############################################################################
# Tests for index_and_collect_contigs

_TEST_FASTA_1_A = FASTA("seq1", "meta1", "TCTTTCAGTCTGGAGACTAGCCTCC")
_TEST_FASTA_1_B = FASTA("seq1", None, "ATTGAGGCGTATTGTGTCG")
_TEST_FASTA_2 = FASTA("seq2", None, "CAAAGCA")
_TEST_FASTA_3 = FASTA("seq3", "more meta data", "AGCTCTCCTCCCC")


def test_index_and_collect_contigs(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    with fasta_file.open("wt") as handle:
        _TEST_FASTA_1_A.write(handle)
        _TEST_FASTA_2.write(handle)
        _TEST_FASTA_3.write(handle)

    assert FASTA.index_and_collect_contigs(fasta_file) == {
        "seq1": 25,
        "seq2": 7,
        "seq3": 13,
    }


def test_index_and_collect_contigs__duplicate_names(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    with fasta_file.open("wt") as handle:
        _TEST_FASTA_1_A.write(handle)
        _TEST_FASTA_2.write(handle)
        _TEST_FASTA_1_B.write(handle)

    assert FASTA.index_and_collect_contigs(fasta_file) == {
        "seq1": 25,
        "seq2": 7,
    }


def test_index_and_collect_contigs__fai_files(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    with fasta_file.open("wt") as handle:
        _TEST_FASTA_1_A.write(handle)

    fai_file = tmp_path / "test.fasta.fai"

    # Fai file should be created once, and then not modified
    FASTA.index_and_collect_contigs(fasta_file)
    stats_1 = fai_file.stat()
    FASTA.index_and_collect_contigs(fasta_file)
    stats_2 = fai_file.stat()

    assert stats_1 == stats_2
