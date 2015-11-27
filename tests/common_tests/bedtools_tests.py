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
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Ignore "invalid" function names
# pylint: disable=C0103
# No need for docstrings
# pylint: disable=C0111
import copy

from nose.tools import \
    assert_equal, \
    assert_not_equal, \
    assert_raises

from paleomix.common.bedtools import \
    BEDRecord


###############################################################################
###############################################################################
# BEDRecord constructor

def test_bedrecord__constructor__defaults():
    record = BEDRecord()
    assert_equal(len(record), 0)
    assert_equal(str(record), "")
    assert_equal(repr(record), "BEDRecord()")


def test_bedrecord__constructor__empty_string():
    record = BEDRecord("")
    assert_equal(len(record), 0)
    assert_equal(str(record), "")
    assert_equal(repr(record), "BEDRecord()")


def test_bedrecord__constructor__3_fields():
    text = "my_contig\t12\t345"
    record = BEDRecord(text)

    assert_equal(len(record), 3)
    assert_equal(str(record), text)
    assert_equal(repr(record),
                 "BEDRecord(contig='my_contig', start=12, end=345)")


def test_bedrecord__constructor__6_fields():
    text = "my_contig\t12\t345\tmy_name\t-3\t-"
    record = BEDRecord(text)

    assert_equal(len(record), 6)
    assert_equal(str(record), text)
    assert_equal(repr(record),
                 "BEDRecord(contig='my_contig', start=12, "
                 "end=345, name='my_name', score=-3, strand='-')")


def test_bedrecord__constructor__extra_fields():
    text = "my_contig\t12\t345\tmy_name\t-3\t-\tfoo\tbar"
    record = BEDRecord(text)

    assert_equal(len(record), 8)
    assert_equal(str(record), text)
    assert_equal(repr(record),
                 "BEDRecord(contig='my_contig', start=12, "
                 "end=345, name='my_name', score=-3, strand='-', "
                 "'foo', 'bar')")


###############################################################################
###############################################################################
# BEDRecord accessors

def test_bedrecord__accessors__3_fields():
    record = BEDRecord("my_contig\t12\t345")

    assert_equal(record.contig, "my_contig")
    assert_equal(record.start, 12)
    assert_equal(record.end, 345)
    assert_raises(IndexError, lambda: record.name)
    assert_raises(IndexError, lambda: record.score)
    assert_raises(IndexError, lambda: record.strand)


def test_bedrecord__accessors__6_fields():
    record = BEDRecord("my_contig\t12\t345\tmy_name\t-3\t-")

    assert_equal(record.contig, "my_contig")
    assert_equal(record.start, 12)
    assert_equal(record.end, 345)
    assert_equal(record.name, "my_name")
    assert_equal(record.score, -3)
    assert_equal(record.strand, "-")


def test_bedrecord__accessors__extra_fields():
    text = "my_contig\t12\t345\tmy_name\t-3\t-\tfoo\tbar"
    record = BEDRecord(text)

    assert_equal(record[6], "foo")
    assert_equal(record[7], "bar")


def test_bedrecord__accessors__6_fields__getitem():
    record = BEDRecord("my_contig\t12\t345\tmy_name\t-3\t-")

    assert_equal(record[0], "my_contig")
    assert_equal(record[1], 12)
    assert_equal(record[2], 345)
    assert_equal(record[3], "my_name")
    assert_equal(record[4], -3)
    assert_equal(record[5], "-")


def test_bedrecord__setters__3_fields():
    record = BEDRecord("my_contig\t12\t345")

    record.contig = "chrZ"
    assert_equal(record.contig, "chrZ")

    record.end += 20
    assert_equal(record.end, 365)

    assert_equal(str(record), "chrZ\t12\t365")
    assert_equal(repr(record),
                 "BEDRecord(contig='chrZ', start=12, end=365)")


def test_bedrecord__setters__type_errors():
    record = BEDRecord("my_contig\t12\t345\tname\t0\t+")

    assert_raises(ValueError, lambda: setattr(record, "contig", 17))
    assert_raises(ValueError, lambda: setattr(record, "start", "foo"))
    assert_raises(ValueError, lambda: setattr(record, "end", "foo"))
    assert_raises(ValueError, lambda: setattr(record, "name", 17.3))
    assert_raises(ValueError, lambda: setattr(record, "score", "foo"))
    assert_raises(ValueError, lambda: setattr(record, "strand", "foo"))


def test_bedrecord__setters__unset_fields__at_end():
    record = BEDRecord("my_contig\t12\t345")

    record.name = "my_region"
    assert_equal(record.name, "my_region")

    record.score = -13
    assert_equal(record.score, -13)

    record.strand = '-'
    assert_equal(record.strand, '-')

    assert_equal(str(record), "my_contig\t12\t345\tmy_region\t-13\t-")
    assert_equal(repr(record),
                 "BEDRecord(contig='my_contig', start=12, end=345, "
                 "name='my_region', score=-13, strand='-')")


def test_bedrecord__setters__unset_fields__after_end():
    record = BEDRecord("")
    record.strand = "-"
    assert_equal(str(record), "\t0\t0\t\t0\t-")

    record = BEDRecord("my_name")
    record.strand = "-"
    assert_equal(str(record), "my_name\t0\t0\t\t0\t-")

    record = BEDRecord("my_name\t17")
    record.strand = "-"
    assert_equal(str(record), "my_name\t17\t0\t\t0\t-")

    record = BEDRecord("my_name\t17\t258")
    record.strand = "-"
    assert_equal(str(record), "my_name\t17\t258\t\t0\t-")

    record = BEDRecord("my_name\t17\t258\tregion")
    record.strand = "-"
    assert_equal(str(record), "my_name\t17\t258\tregion\t0\t-")

    record = BEDRecord("my_name\t17\t258\tregion\t33")
    record.strand = "-"
    assert_equal(str(record), "my_name\t17\t258\tregion\t33\t-")

    record = BEDRecord("my_name\t17\t258\tregion\t33\t+")
    record.strand = "-"
    assert_equal(str(record), "my_name\t17\t258\tregion\t33\t-")


def test_bedrecord__cmp():
    record_1_txt = "my_contig\t12\t345\tmy_name\t-3\t-\tfoo"
    record_1 = BEDRecord(record_1_txt)
    record_2 = BEDRecord("chrZ\t132\t4345\tchrZ_region\t0\t+\tbar")

    for idx in xrange(len(record_2)):
        record_tmp = BEDRecord(record_1_txt)
        assert_equal(record_1, record_tmp)
        record_tmp[idx] = record_2[idx]
        assert_not_equal(record_1, record_tmp)
        record_tmp[idx] = record_1[idx]
        assert_equal(record_1, record_tmp)


###############################################################################
###############################################################################

def test_bedrecord__copy():
    record_1_txt = "my_contig\t12\t345\tmy_name\t-3\t-"
    record_1 = BEDRecord(record_1_txt)
    record_2 = copy.copy(record_1)
    record_2.name = "my_clone"

    assert_equal(str(record_1), record_1_txt)
    assert_equal(str(record_2), "my_contig\t12\t345\tmy_clone\t-3\t-")


def test_bedrecord__deepcopy():
    record_1_txt = "my_contig\t12\t345\tmy_name\t-3\t-"
    record_1 = BEDRecord(record_1_txt)
    record_1[6] = ["foo"]
    record_2 = copy.deepcopy(record_1)
    record_2[6][0] = "bar"

    assert_equal(str(record_1), record_1_txt + "\t['foo']")
    assert_equal(str(record_2), record_1_txt + "\t['bar']")
