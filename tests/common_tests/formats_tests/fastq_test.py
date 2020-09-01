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
from paleomix.common.formats.fastq import FASTQ, FASTQError, FASTQualities


_SEQ_FRAG = "AAGTCC"  # len() = 6
_QUAL_FRAG = "123456"  # len() = 6


###############################################################################
###############################################################################
# Tests for FASTQ constructor


def _simple_fastq_record():
    return FASTQ("Dummy", "Meta-inf", "ACGT", "1234")


def test_fastq__constructor():
    record = _simple_fastq_record()
    assert record.name == "Dummy"
    assert record.meta == "Meta-inf"
    assert record.sequence == "ACGT"
    assert record.qualities == "1234"


def test_fastq__constructor__no_meta():
    record = FASTQ("Dummy", None, "ACGT", "1234")
    assert record.name == "Dummy"
    assert record.meta == ""
    assert record.sequence == "ACGT"
    assert record.qualities == "1234"


def test_fastq__constructor__name_must_be_string_type():
    with pytest.raises(FASTQError, match="FASTQ name must be a non-empty string"):
        FASTQ(1, None, "ACGT", "1234")


def test_fastq__constructor__name_must_be_non_empty():
    with pytest.raises(FASTQError, match="FASTQ name must be a non-empty string"):
        FASTQ("", None, "ACGT", "1234")


def test_fastq__constructor__meta_must_be_string_type_or_none():
    with pytest.raises(FASTQError, match="FASTQ meta must be a string, or None"):
        FASTQ("Seq1", 1, "ACGT", "1234")


def test_fastq__constructor__seq_must_be_string_type():
    with pytest.raises(FASTQError, match="FASTQ sequence must be a string"):
        FASTQ("Seq1", None, 1, "1234")


def test_fastq__constructor__qual_must_be_string_type():
    with pytest.raises(FASTQError, match="FASTQ qualities must be a string"):
        FASTQ("Seq1", None, "ACGT", 1234)


def test_fastq__constructor__seq_and_qual_must_have_same_len():
    with pytest.raises(FASTQError, match="FASTQ sequence and qualities length differ"):
        FASTQ("Seq1", None, "ACGT", "12345")


###############################################################################
###############################################################################
# Tests for __repr__


def test_fastq__repr():
    expected = "FASTQ('foobar', '', %r, %r)" % (_SEQ_FRAG, _QUAL_FRAG)
    result = repr(FASTQ("foobar", None, _SEQ_FRAG, _QUAL_FRAG))
    assert result == expected


def test_fastq__repr_with_meta_information():
    expected = "FASTQ('foobar', 'my Meta-Info', %r, %r)" % (_SEQ_FRAG, _QUAL_FRAG)
    result = repr(FASTQ("foobar", "my Meta-Info", _SEQ_FRAG, _QUAL_FRAG))
    assert result == expected


###############################################################################
###############################################################################
# Tests for write


def test_fastq__write():
    expected = "@foobar\n%s\n+\n%s\n" % (_SEQ_FRAG, _QUAL_FRAG)
    stringf = io.StringIO()
    FASTQ("foobar", None, _SEQ_FRAG, _QUAL_FRAG).write(stringf)
    assert stringf.getvalue() == expected


def test_fastq__write_with_meta_information():
    expected = "@foobar my Meta-Info\n%s\n+\n%s\n" % (_SEQ_FRAG, _QUAL_FRAG)
    stringf = io.StringIO()
    FASTQ("foobar", "my Meta-Info", _SEQ_FRAG, _QUAL_FRAG).write(stringf)
    assert stringf.getvalue() == expected


###############################################################################
###############################################################################
# Tests for FASTQ.from_lines


def test_fasta__from_lines__no_records():
    assert list(FASTQ.from_lines([])) == list([])


def test_fasta__from_lines__no_records__empty_line():
    assert list(FASTQ.from_lines([""])) == list([])


def test_fasta__from_lines_single_record():
    lines = ["@single\n", "CCTTCATCCA\n", "+", "1234567890"]
    expected = [FASTQ("single", None, "CCTTCATCCA", "1234567890")]
    assert list(FASTQ.from_lines(lines)) == list(expected)


def test_fasta__from_lines__multiple_records():
    lines = [
        "@first\n",
        "TGTTCTCCACCGTGCACAAC\n",
        "+",
        "12345678901234567890\n",
        "@Second XT:1:0\n",
        "GAGAGCTCAGCTAAC\n",
        "+\n",
        "098765432109876\n",
        "@Third\n",
        "GGCATTCGGC\n",
        "+\n",
        "5678901234\n",
    ]
    expected = [
        FASTQ("first", None, "TGTTCTCCACCGTGCACAAC", "12345678901234567890"),
        FASTQ("Second", "XT:1:0", "GAGAGCTCAGCTAAC", "098765432109876"),
        FASTQ("Third", None, "GGCATTCGGC", "5678901234"),
    ]
    assert list(FASTQ.from_lines(lines)) == list(expected)


def test_fasta__from_lines__partial__1_line():
    with pytest.raises(FASTQError, match="Partial FASTQ record"):
        list(FASTQ.from_lines(["@fastq1\n"]))


def test_fasta__from_lines__partial__2_lines():
    with pytest.raises(FASTQError, match="Partial FASTQ record"):
        list(FASTQ.from_lines(["@fastq1\n", "ACGT\n"]))


def test_fasta__from_lines__partial__3_lines():
    with pytest.raises(FASTQError, match="Partial FASTQ record"):
        list(FASTQ.from_lines(["@fastq1\n", "ACGT\n", "+\n"]))


def test_fasta__from_lines__invalid_header():
    lines = [">fastq\n", "GGCATTCGGC\n", "+\n", "5678901234\n"]
    with pytest.raises(FASTQError, match="Invalid FASTQ header"):
        list(FASTQ.from_lines(lines))


def test_fasta__from_lines__mismatching_lengths():
    lines = ["@fastq\n", "GGCATTCGGC\n", "+\n", "567890123\n"]
    with pytest.raises(
        FASTQError, match="Sequence length does not match qualities length"
    ):
        list(FASTQ.from_lines(lines))


def test_fasta__from_lines__invalid_separator():
    lines = ["@fastq\n", "GGCATTCGGC\n", "?\n", "5678901234\n"]
    with pytest.raises(FASTQError, match="Invalid FASTQ separator"):
        list(FASTQ.from_lines(lines))


###############################################################################
###############################################################################
# Tests for 'FASTQ.from_file'


@pytest.mark.parametrize("func", (open, gzip.open, bz2.open))
def test_fasta__from_file(func, tmp_path):
    expected = [
        FASTQ("This_is_FASTA!", None, "ACGTN", "12345"),
        FASTQ("This_is_ALSO_FASTA!", None, "CGTNA", "56789"),
    ]

    with func(fspath(tmp_path / "file"), "wt") as handle:
        for item in expected:
            item.write(handle)

    assert list(FASTQ.from_file(tmp_path / "file")) == expected


###############################################################################
###############################################################################
# Tests for 'FASTQualities'

_33_READ = FASTQ("33", None, "ACGT", "!02I")
_64_READ = FASTQ("33", None, "ACGT", "@UXi")
_AMBIGIOUS_read = FASTQ("33", None, "ACGT", ";CDI")


def test_fastqualities__no_qualities():
    quals = FASTQualities()

    assert quals.offsets() == FASTQualities.MISSING


def test_fastqualities__phred_33():
    quals = FASTQualities()
    quals.update(_33_READ)

    assert quals.offsets() == FASTQualities.OFFSET_33


def test_fastqualities__phred_64():
    quals = FASTQualities()
    quals.update(_64_READ)

    assert quals.offsets() == FASTQualities.OFFSET_64


def test_fastqualities__phred_ambigious():
    quals = FASTQualities()
    quals.update(_AMBIGIOUS_read)

    assert quals.offsets() == FASTQualities.AMBIGIOUS


def test_fastqualities__phred_both():
    quals = FASTQualities()
    quals.update(_33_READ)
    quals.update(_64_READ)

    assert quals.offsets() == FASTQualities.BOTH


###############################################################################
###############################################################################


def test_fasta__equality():
    assert FASTQ("A", "B", "C", "D") == FASTQ("A", "B", "C", "D")


def test_fasta__inequality():
    assert FASTQ("A", "B", "C", "D") != FASTQ("A", "B", "D", "D")
    assert FASTQ("A", "B", "C", "D") != FASTQ("A", None, "C", "D")
    assert FASTQ("A", "B", "C", "D") != FASTQ("D", "B", "C", "D")


def test_fasta__sorting_less_equal():
    assert not FASTQ("A", "B", "C", "D") < FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "B", "C", "D") < FASTQ("B", "B", "C", "D")
    assert FASTQ("A", "B", "C", "D") < FASTQ("A", "C", "C", "D")
    assert FASTQ("A", "B", "C", "D") < FASTQ("A", "B", "D", "D")
    assert FASTQ("A", "B", "C", "D") <= FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "B", "C", "D") <= FASTQ("B", "B", "C", "D")
    assert FASTQ("A", "B", "C", "D") <= FASTQ("A", "C", "C", "D")
    assert FASTQ("A", "B", "C", "D") <= FASTQ("A", "B", "D", "D")


def test_fasta__sorting_greater_equal():
    assert not FASTQ("A", "B", "C", "D") > FASTQ("A", "B", "C", "D")
    assert FASTQ("B", "B", "C", "D") > FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "C", "C", "D") > FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "B", "D", "D") > FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "B", "C", "E") > FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "B", "C", "D") >= FASTQ("A", "B", "C", "D")
    assert FASTQ("B", "B", "C", "D") >= FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "C", "C", "D") >= FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "B", "D", "D") >= FASTQ("A", "B", "C", "D")
    assert FASTQ("A", "B", "C", "E") >= FASTQ("A", "B", "C", "D")


def test_fasta__hash():
    assert hash(FASTQ("A", "B", "C", "D")) == hash(FASTQ("A", "B", "C", "D"))
    assert hash(FASTQ("A", "B", "C", "D")) != hash(FASTQ("B", "B", "C", "D"))
    assert hash(FASTQ("A", "B", "C", "D")) != hash(FASTQ("A", "C", "C", "D"))
    assert hash(FASTQ("A", "B", "C", "D")) != hash(FASTQ("A", "B", "D", "D"))
    assert hash(FASTQ("A", "B", "C", "D")) != hash(FASTQ("A", "B", "C", "E"))


def test_fasta__unimplemented_comparison():
    assert NotImplemented is FASTQ("A", None, "C", "D").__eq__(10)
    assert NotImplemented is FASTQ("A", None, "C", "D").__lt__(10)
    assert NotImplemented is FASTQ("A", None, "C", "D").__le__(10)
    assert NotImplemented is FASTQ("A", None, "C", "D").__ge__(10)
    assert NotImplemented is FASTQ("A", None, "C", "D").__gt__(10)
