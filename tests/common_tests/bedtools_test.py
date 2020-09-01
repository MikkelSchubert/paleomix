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
import copy

import pytest

from paleomix.common.bedtools import (
    BEDError,
    BEDRecord,
    merge_bed_records,
    pad_bed_records,
)

###############################################################################
###############################################################################
# BEDRecord constructor


def test_bedrecord__constructor__defaults():
    record = BEDRecord()
    assert len(record) == 0
    assert str(record) == ""
    assert repr(record) == "BEDRecord()"


def test_bedrecord__constructor__empty_string():
    record = BEDRecord("")
    assert len(record) == 0
    assert str(record) == ""
    assert repr(record) == "BEDRecord()"


def test_bedrecord__constructor__3_fields():
    text = "my_contig\t12\t345"
    record = BEDRecord(text)

    assert len(record) == 3
    assert str(record) == text
    assert repr(record) == "BEDRecord(contig='my_contig', start=12, end=345)"


def test_bedrecord__constructor__6_fields():
    text = "my_contig\t12\t345\tmy_name\t-3\t-"
    record = BEDRecord(text)

    assert len(record) == 6
    assert str(record) == text
    assert (
        repr(record) == "BEDRecord(contig='my_contig', start=12, "
        "end=345, name='my_name', score=-3, strand='-')"
    )


def test_bedrecord__constructor__extra_fields():
    text = "my_contig\t12\t345\tmy_name\t-3\t-\tfoo\tbar"
    record = BEDRecord(text)

    assert len(record) == 8
    assert str(record) == text
    assert (
        repr(record) == "BEDRecord(contig='my_contig', start=12, "
        "end=345, name='my_name', score=-3, strand='-', "
        "'foo', 'bar')"
    )


def test_bedrecord__constructor__invalid_values():
    tmpl = "\t1\t%s\t\t0\t-"
    BEDRecord(tmpl % (1,))  # check template

    with pytest.raises(BEDError):
        BEDRecord(tmpl % ("not a number",))


###############################################################################
###############################################################################
# BEDRecord accessors


def test_bedrecord__accessors__3_fields():
    record = BEDRecord("my_contig\t12\t345")

    assert record.contig == "my_contig"
    assert record.start == 12
    assert record.end == 345
    with pytest.raises(IndexError):
        record.name()
    with pytest.raises(IndexError):
        record.score()
    with pytest.raises(IndexError):
        record.strand()


def test_bedrecord__accessors__6_fields():
    record = BEDRecord("my_contig\t12\t345\tmy_name\t-3\t-")

    assert record.contig == "my_contig"
    assert record.start == 12
    assert record.end == 345
    assert record.name == "my_name"
    assert record.score == -3
    assert record.strand == "-"


def test_bedrecord__accessors__extra_fields():
    text = "my_contig\t12\t345\tmy_name\t-3\t-\tfoo\tbar"
    record = BEDRecord(text)

    assert record[6] == "foo"
    assert record[7] == "bar"


def test_bedrecord__accessors__6_fields__getitem():
    record = BEDRecord("my_contig\t12\t345\tmy_name\t-3\t-")

    assert record[0] == "my_contig"
    assert record[1] == 12
    assert record[2] == 345
    assert record[3] == "my_name"
    assert record[4] == -3
    assert record[5] == "-"


def test_bedrecord__setters__3_fields():
    record = BEDRecord("my_contig\t12\t345")

    record.contig = "chrZ"
    assert record.contig == "chrZ"

    record.end += 20
    assert record.end == 365

    assert str(record) == "chrZ\t12\t365"
    assert repr(record) == "BEDRecord(contig='chrZ', start=12, end=365)"


def test_bedrecord__setters__type_errors():
    record = BEDRecord("my_contig\t12\t345\tname\t0\t+")

    with pytest.raises(ValueError):
        record.contig = 17
    with pytest.raises(ValueError):
        record.start = "foo"
    with pytest.raises(ValueError):
        record.end = "foo"
    with pytest.raises(ValueError):
        record.name = 17.3
    with pytest.raises(ValueError):
        record.score = "foo"
    with pytest.raises(ValueError):
        record.strand = "foo"


def test_bedrecord__setters__unset_fields__at_end():
    record = BEDRecord("my_contig\t12\t345")

    record.name = "my_region"
    assert record.name == "my_region"

    record.score = -13
    assert record.score == -13

    record.strand = "-"
    assert record.strand == "-"

    assert str(record) == "my_contig\t12\t345\tmy_region\t-13\t-"
    assert (
        repr(record) == "BEDRecord(contig='my_contig', start=12, end=345, "
        "name='my_region', score=-13, strand='-')"
    )


def test_bedrecord__setters__unset_fields__after_end():
    record = BEDRecord("")
    record.strand = "-"
    assert str(record) == "\t0\t0\t\t0\t-"

    record = BEDRecord("my_name")
    record.strand = "-"
    assert str(record) == "my_name\t0\t0\t\t0\t-"

    record = BEDRecord("my_name\t17")
    record.strand = "-"
    assert str(record) == "my_name\t17\t0\t\t0\t-"

    record = BEDRecord("my_name\t17\t258")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\t\t0\t-"

    record = BEDRecord("my_name\t17\t258\tregion")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\tregion\t0\t-"

    record = BEDRecord("my_name\t17\t258\tregion\t33")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\tregion\t33\t-"

    record = BEDRecord("my_name\t17\t258\tregion\t33\t+")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\tregion\t33\t-"


def test_bedrecord__cmp():
    record_1_txt = "my_contig\t12\t345\tmy_name\t-3\t-\tfoo"
    record_1 = BEDRecord(record_1_txt)
    record_2 = BEDRecord("chrZ\t132\t4345\tchrZ_region\t0\t+\tbar")

    for idx in range(len(record_2)):
        record_tmp = BEDRecord(record_1_txt)
        assert record_1 == record_tmp
        record_tmp[idx] = record_2[idx]
        assert record_1 != record_tmp
        record_tmp[idx] = record_1[idx]
        assert record_1 == record_tmp


###############################################################################
# freeze


def test_bedrecord__freeze__modify_existing_value():
    record_str = "my_contig\t12\t345\tmy_name"
    record = BEDRecord(record_str)
    frozen_record = record.freeze()

    with pytest.raises(TypeError):
        frozen_record.start = 1

    assert str(frozen_record) == record_str


def test_bedrecord__freeze__add_new_value():
    record_str = "my_contig\t12\t345\tmy_name"
    record = BEDRecord(record_str)
    frozen_record = record.freeze()

    with pytest.raises(AttributeError):
        frozen_record.score = 1

    assert str(frozen_record) == record_str


###############################################################################
###############################################################################


def test_bedrecord__copy():
    record_1_txt = "my_contig\t12\t345\tmy_name\t-3\t-"
    record_1 = BEDRecord(record_1_txt)
    record_2 = copy.copy(record_1)
    record_2.name = "my_clone"

    assert str(record_1) == record_1_txt
    assert str(record_2) == "my_contig\t12\t345\tmy_clone\t-3\t-"


def test_bedrecord__deepcopy():
    record_1_txt = "my_contig\t12\t345\tmy_name\t-3\t-"
    record_1 = BEDRecord(record_1_txt)
    record_1[6] = ["foo"]
    record_2 = copy.deepcopy(record_1)
    record_2[6][0] = "bar"

    assert str(record_1) == record_1_txt + "\t['foo']"
    assert str(record_2) == record_1_txt + "\t['bar']"


###############################################################################
# pad_bed_records


def test_pad_bed_records__empty_sequences():
    pad_bed_records([], 10)
    pad_bed_records((), 10)


def test_pad_bed_records():
    records = [
        _new_bed_record("chr1", 10, 90),
        _new_bed_record("chr2", 100, 200),
    ]

    pad_bed_records(records, 20)
    assert records == [
        _new_bed_record("chr1", 0, 110),
        _new_bed_record("chr2", 80, 220),
    ]


def test_pad_bed_records__with_max_lengths():
    max_sizes = {"chr1": 100, "chr2": 200}
    records = [
        _new_bed_record("chr1", 10, 90),
        _new_bed_record("chr2", 10, 90),
        _new_bed_record("chr2", 100, 190),
    ]

    pad_bed_records(records, 20, max_sizes)
    assert records == [
        _new_bed_record("chr1", 0, 100),
        _new_bed_record("chr2", 0, 110),
        _new_bed_record("chr2", 80, 200),
    ]


###############################################################################
# merge_bed_records


def test_merge_records__empty_sequences():
    assert merge_bed_records(()) == []
    assert merge_bed_records([]) == []


def test_merge_records__single_record():
    assert merge_bed_records([_new_bed_record("chr1", 1234, 5678)]) == [
        _new_bed_record("chr1", 1234, 5678)
    ]


def test_merge_records__minimal_fields_only():
    assert merge_bed_records([_new_bed_record("chr1", 1234, 5678, "foo", 1, "-")]) == [
        _new_bed_record("chr1", 1234, 5678)
    ]


def test_merge_records__overlapping_records_1():
    assert merge_bed_records(
        [_new_bed_record("chr1", 1234, 5678), _new_bed_record("chr1", 5677, 9012)]
    ) == [_new_bed_record("chr1", 1234, 9012)]


def test_merge_records__overlapping_records_2():
    assert merge_bed_records(
        [_new_bed_record("chr1", 1234, 5678), _new_bed_record("chr1", 5678, 9012)]
    ) == [_new_bed_record("chr1", 1234, 9012)]


def test_merge_records__non_overlapping_records_1():
    assert merge_bed_records(
        [_new_bed_record("chr1", 1234, 5678), _new_bed_record("chr1", 5679, 9012)]
    ) == [_new_bed_record("chr1", 1234, 5678), _new_bed_record("chr1", 5679, 9012)]


def test_merge_records__non_overlapping_records_2():
    assert merge_bed_records(
        [_new_bed_record("chr1", 1234, 5678), _new_bed_record("chr1", 5680, 9012)]
    ) == [_new_bed_record("chr1", 1234, 5678), _new_bed_record("chr1", 5680, 9012)]


def test_merge_records__complex_example():
    assert merge_bed_records(
        [
            _new_bed_record("chr1", 1234, 5678),
            _new_bed_record("chr1", 5678, 9012),
            _new_bed_record("chr2", 1, 20),
            _new_bed_record("chr2", 100, 200),
            _new_bed_record("chr2", 150, 250),
        ]
    ) == [
        _new_bed_record("chr1", 1234, 9012),
        _new_bed_record("chr2", 1, 20),
        _new_bed_record("chr2", 100, 250),
    ]


def test_merge_records__complex_example__unsorted():
    assert merge_bed_records(
        [
            _new_bed_record("chr2", 100, 200),
            _new_bed_record("chr1", 1234, 5678),
            _new_bed_record("chr2", 150, 250),
            _new_bed_record("chr2", 1, 20),
            _new_bed_record("chr1", 5678, 9012),
        ]
    ) == [
        _new_bed_record("chr1", 1234, 9012),
        _new_bed_record("chr2", 1, 20),
        _new_bed_record("chr2", 100, 250),
    ]


def _new_bed_record(*args):
    record = BEDRecord()
    record._fields = list(args)

    return record
