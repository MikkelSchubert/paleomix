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
from pathlib import Path
from typing import Any, List
from unittest.mock import Mock

import pytest

from paleomix.common.formats.bed import (
    BEDError,
    BEDRecord,
    merge_bed_records,
    pad_bed_records,
    read_bed_file,
    sort_bed_by_bamfile,
)

########################################################################################
# BEDRecord constructor


def test_bedrecord__constructor_3():
    record = BEDRecord("my_contig", 12, 345)
    assert record.contig == "my_contig"
    assert record.start == 12
    assert record.end == 345
    assert record.name is None
    assert record.score is None
    assert record.strand is None
    assert str(record) == "my_contig\t12\t345"
    assert repr(record) == "BEDRecord(contig='my_contig', start=12, end=345)"


def test_bedrecord__constructor_6():
    record = BEDRecord("my_contig", 12, 345, "my_name", -3, "-")
    assert record.contig == "my_contig"
    assert record.start == 12
    assert record.end == 345
    assert record.name == "my_name"
    assert record.score == -3
    assert record.strand == "-"
    assert str(record) == "my_contig\t12\t345\tmy_name\t-3\t-"
    assert (
        repr(record) == "BEDRecord(contig='my_contig', start=12, "
        "end=345, name='my_name', score=-3, strand='-')"
    )


@pytest.mark.parametrize("contig", (None, ""))
def test_bedrecord__constructor__empty_contig(contig: Any):
    with pytest.raises(ValueError, match="contig is blank"):
        BEDRecord(contig, 12, 345, "my_name", -3, "-")


@pytest.mark.parametrize("strand", ("", "?", "foo"))
def test_bedrecord__constructor__invalid_strand(strand: Any):
    with pytest.raises(ValueError, match="invalid strand"):
        BEDRecord("contig", 12, 345, "my_name", -3, strand)


########################################################################################
# BEDRecord.parse


@pytest.mark.parametrize("text", ("", "contig", "contig\t0"))
def test_bedrecord__parse__0_to_2_fields(text: str):
    with pytest.raises(BEDError, match="not enough columns"):
        BEDRecord.parse(text)


def test_bedrecord__parse__3_fields():
    text = "my_contig\t12\t345"
    record = BEDRecord.parse(text)

    assert str(record) == text
    assert repr(record) == "BEDRecord(contig='my_contig', start=12, end=345)"


_6_COLUMN_LINES = [
    "my_contig\t12\t345\tmy_name\t-3\t-",
    "my_contig\t12\t345\tmy_name\t-3\t-\tfoo\tbar",
]


@pytest.mark.parametrize("text", _6_COLUMN_LINES)
def test_bedrecord__parse__6_fields(text: str):
    record = BEDRecord.parse(text)
    assert record.contig == "my_contig"
    assert record.start == 12
    assert record.end == 345
    assert record.name == "my_name"
    assert record.score == -3
    assert record.strand == "-"
    assert str(record) == "\t".join(text.split("\t")[:6])
    assert (
        repr(record) == "BEDRecord(contig='my_contig', start=12, "
        "end=345, name='my_name', score=-3, strand='-')"
    )


def test_bedrecord__parse__invalid_values_1():
    with pytest.raises(BEDError, match="contig is blank"):
        BEDRecord.parse("\t123\t456")


def test_bedrecord__parse__invalid_values_2():
    with pytest.raises(BEDError, match="Expected int in column 2"):
        BEDRecord.parse("contig\tsix\t456")


def test_bedrecord__parse__invalid_values_3():
    with pytest.raises(BEDError, match="Expected int in column 3"):
        BEDRecord.parse("contig\t123\tlots")


def test_bedrecord__parse__invalid_values_5():
    with pytest.raises(BEDError, match="Expected int in column 5"):
        BEDRecord.parse("contig\t123\t456\tfoo\tbar")


def test_bedrecord__parse__invalid_values_6():
    with pytest.raises(BEDError, match="strand must be \\+ or -"):
        BEDRecord.parse("contig\t123\t456\tfoo\t0\t?")


def test_bedrecord__parse__invalid_values_1_():
    tmpl = "contig\t0\t%s\t\t0\t-"
    BEDRecord.parse(tmpl % (1,))  # check template

    with pytest.raises(BEDError, match=""):
        BEDRecord.parse(tmpl % ("not a number",))


def test_bedrecord__setters__unset_fields__after_end():
    record = BEDRecord.parse("my_name\t17\t258")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\t\t0\t-"

    record = BEDRecord.parse("my_name\t17\t258\tregion")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\tregion\t0\t-"

    record = BEDRecord.parse("my_name\t17\t258\tregion\t33")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\tregion\t33\t-"

    record = BEDRecord.parse("my_name\t17\t258\tregion\t33\t+")
    record.strand = "-"
    assert str(record) == "my_name\t17\t258\tregion\t33\t-"


###############################################################################
# BEDRecord comparisons


def test_sorting__compared_with_non_bed_record():
    with pytest.raises(TypeError):
        _ = BEDRecord("chr2", 1, 20) < "foo"


def test_bedrecord__cmp():
    record_1 = BEDRecord("my_contig", 12, 345, "my_name", -3, "-")
    record_2 = BEDRecord("other", 565, 684, "name2", 0, "+")

    for key in ("contig", "start", "end", "name", "score", "strand"):
        record_tmp = copy.copy(record_1)
        assert record_1 == record_tmp
        setattr(record_tmp, key, getattr(record_2, key))
        assert record_1 != record_tmp
        setattr(record_tmp, key, getattr(record_1, key))
        assert record_1 == record_tmp


###############################################################################
# pad_bed_records


def test_pad_bed_records__empty_sequences():
    pad_bed_records([], 10)
    pad_bed_records((), 10)


def test_pad_bed_records():
    records = [
        BEDRecord("chr1", 10, 90),
        BEDRecord("chr2", 100, 200),
    ]

    assert pad_bed_records(records, 20) == [
        BEDRecord("chr1", 0, 110),
        BEDRecord("chr2", 80, 220),
    ]


def test_pad_bed_records__with_max_lengths():
    max_sizes = {"chr1": 100, "chr2": 200}
    records = [
        BEDRecord("chr1", 10, 90),
        BEDRecord("chr2", 10, 90),
        BEDRecord("chr2", 100, 190),
    ]

    assert pad_bed_records(records, 20, max_sizes) == [
        BEDRecord("chr1", 0, 100),
        BEDRecord("chr2", 0, 110),
        BEDRecord("chr2", 80, 200),
    ]


def test_pad_bed_records__negative_padding():
    records = [
        BEDRecord("chr1", 10, 90),
        BEDRecord("chr2", 100, 200),
    ]

    assert pad_bed_records(records, -15) == [
        BEDRecord("chr1", 25, 75),
        BEDRecord("chr2", 115, 185),
    ]


def test_pad_bed_records__negative_padding__near_empty_records():
    assert pad_bed_records([BEDRecord("chr1", 10, 90)], -38) == [
        BEDRecord("chr1", 48, 52)
    ]
    assert pad_bed_records([BEDRecord("chr1", 10, 90)], -39) == [
        BEDRecord("chr1", 49, 51)
    ]
    assert pad_bed_records([BEDRecord("chr1", 10, 91)], -40) == [
        BEDRecord("chr1", 50, 51)
    ]


def test_pad_bed_records__negative_padding__empty_records():
    assert pad_bed_records([BEDRecord("chr1", 10, 90)], -40) == []
    assert pad_bed_records([BEDRecord("chr1", 10, 90)], -41) == []
    assert pad_bed_records([BEDRecord("chr1", 10, 91)], -41) == []


###############################################################################
# merge_bed_records


def test_merge_records__empty_sequences():
    assert merge_bed_records(()) == []
    assert merge_bed_records([]) == []


def test_merge_records__single_record():
    assert merge_bed_records([BEDRecord("chr1", 1234, 5678)]) == [
        BEDRecord("chr1", 1234, 5678)
    ]


def test_merge_records__minimal_fields_only():
    assert merge_bed_records([BEDRecord("chr1", 1234, 5678, "foo", 1, "-")]) == [
        BEDRecord("chr1", 1234, 5678)
    ]


def test_merge_records__overlapping_records_1():
    assert merge_bed_records(
        [BEDRecord("chr1", 1234, 5678), BEDRecord("chr1", 5677, 9012)]
    ) == [BEDRecord("chr1", 1234, 9012)]


def test_merge_records__overlapping_records_2():
    assert merge_bed_records(
        [BEDRecord("chr1", 1234, 5678), BEDRecord("chr1", 5678, 9012)]
    ) == [BEDRecord("chr1", 1234, 9012)]


def test_merge_records__non_overlapping_records_1():
    assert merge_bed_records(
        [BEDRecord("chr1", 1234, 5678), BEDRecord("chr1", 5679, 9012)]
    ) == [BEDRecord("chr1", 1234, 5678), BEDRecord("chr1", 5679, 9012)]


def test_merge_records__non_overlapping_records_2():
    assert merge_bed_records(
        [BEDRecord("chr1", 1234, 5678), BEDRecord("chr1", 5680, 9012)]
    ) == [BEDRecord("chr1", 1234, 5678), BEDRecord("chr1", 5680, 9012)]


def test_merge_records__complex_example():
    assert merge_bed_records(
        [
            BEDRecord("chr1", 1234, 5678),
            BEDRecord("chr1", 5678, 9012),
            BEDRecord("chr2", 1, 20),
            BEDRecord("chr2", 100, 200),
            BEDRecord("chr2", 150, 250),
        ]
    ) == [
        BEDRecord("chr1", 1234, 9012),
        BEDRecord("chr2", 1, 20),
        BEDRecord("chr2", 100, 250),
    ]


def test_merge_records__complex_example__unsorted():
    assert merge_bed_records(
        [
            BEDRecord("chr2", 100, 200),
            BEDRecord("chr1", 1234, 5678),
            BEDRecord("chr2", 150, 250),
            BEDRecord("chr2", 1, 20),
            BEDRecord("chr1", 5678, 9012),
        ]
    ) == [
        BEDRecord("chr1", 1234, 9012),
        BEDRecord("chr2", 1, 20),
        BEDRecord("chr2", 100, 250),
    ]


###############################################################################
# read_bed_file

_SIMPLE_BED = """chr1\t7\t9
chr3\t123\t45123\tFoo\t0
chr9\t777\t999
"""

_SIMPLE_BED_WITH_SKIPPED_LINES = """
chr1\t7\t9
# chr3\t123\t45123\tFoo\t0
chr9\t777\t999
"""


def _write_bed(tmp_path: Path, data: str) -> str:
    filename = tmp_path / "tmp.bed"
    with filename.open("wt") as handle:
        handle.write(data)

    return str(filename)


def test_read_bed_file__empty_file(tmp_path: Path):
    filename = _write_bed(tmp_path, "")

    assert list(read_bed_file(filename)) == []


def test_read_bed_file__simple_records(tmp_path: Path):
    filename = _write_bed(tmp_path, _SIMPLE_BED)

    assert list(read_bed_file(filename)) == [
        BEDRecord("chr1", 7, 9),
        BEDRecord("chr3", 123, 45123, "Foo", 0),
        BEDRecord("chr9", 777, 999),
    ]


def test_read_bed_file__simple_records_and_skipped_lines(tmp_path: Path):
    filename = _write_bed(tmp_path, _SIMPLE_BED_WITH_SKIPPED_LINES)

    assert list(read_bed_file(filename)) == [
        BEDRecord("chr1", 7, 9),
        BEDRecord("chr9", 777, 999),
    ]


def test_read_bed_file__parse_error(tmp_path: Path):
    filename = _write_bed(tmp_path, "chr1\t0\tabc")

    with pytest.raises(BEDError, match=":1: Expected int in column 3 but found 'abc'"):
        next(read_bed_file(filename))


_INVALID_START_END = ["chr1\t0\t0", "chr1\t-1\t100", "chr1\t-100\t-1", "chr1\t100\t10"]
_POSITIONS_INSIDE_CONTIG = ["chr1\t0\t50", "chr1\t50\t99", "chr1\t99\t100"]
_POSITIONS_OUTSIDE_CONTIG = ["chr1\t0\t101", "chr1\t99\t102", "chr1\t200\t300"]


@pytest.mark.parametrize("text", _INVALID_START_END)
def test_read_bed_file__invalid_coordinates(tmp_path: Path, text: str):
    filename = _write_bed(tmp_path, text)
    with pytest.raises(BEDError, match="invalid start/end coordinates"):
        list(read_bed_file(filename))


@pytest.mark.parametrize("text", _POSITIONS_INSIDE_CONTIG)
def test_read_bed_file__inside_contig(tmp_path: Path, text: str):
    filename = _write_bed(tmp_path, text)
    list(read_bed_file(filename, {"chr1": 100}))


@pytest.mark.parametrize("text", _POSITIONS_OUTSIDE_CONTIG)
def test_read_bed_file__outside_contig(tmp_path: Path, text: str):
    filename = _write_bed(tmp_path, text)
    with pytest.raises(BEDError, match="coordinates outside contig"):
        list(read_bed_file(filename, {"chr1": 100}))


def test_read_bed_file__unknown_contigs(tmp_path: Path):
    filename = _write_bed(tmp_path, "chr2\t0\t100")
    with pytest.raises(BEDError, match="unknown contig"):
        list(read_bed_file(filename, {"chr1": 100}))


########################################################################################
# sort_bed_by_bamfile


_UNSORTED_RECORDS = [
    BEDRecord("chr3", 100, 200),
    BEDRecord("chr2", 500, 10010),
    BEDRecord("chr3", 1, 10),
    BEDRecord("chr1", 0, 1234),
]


def test_sort_bed_by_bamfile__empty_list():
    handle: Any = Mock()
    handle.references = []
    regions: List[BEDRecord] = []

    sort_bed_by_bamfile(handle, regions)
    assert regions == []


def test_sort_bed_by_bamfile__known_contigs_1():
    handle: Any = Mock()
    handle.references = ["chr1", "chr2", "chr3"]
    regions = list(_UNSORTED_RECORDS)

    sort_bed_by_bamfile(handle, regions)
    assert regions == [
        BEDRecord("chr1", 0, 1234),
        BEDRecord("chr2", 500, 10010),
        BEDRecord("chr3", 1, 10),
        BEDRecord("chr3", 100, 200),
    ]


def test_sort_bed_by_bamfile__known_contigs_2():
    handle: Any = Mock()
    handle.references = ["chr3", "chr2", "chr1"]
    regions = list(_UNSORTED_RECORDS)

    sort_bed_by_bamfile(handle, regions)
    assert regions == [
        BEDRecord("chr3", 1, 10),
        BEDRecord("chr3", 100, 200),
        BEDRecord("chr2", 500, 10010),
        BEDRecord("chr1", 0, 1234),
    ]


def test_sort_bed_by_bamfile__unknown_contigs_1():
    handle: Any = Mock()
    handle.references = ["chr1", "chr2"]
    regions = list(_UNSORTED_RECORDS)

    sort_bed_by_bamfile(handle, regions)
    assert regions == [
        BEDRecord("chr1", 0, 1234),
        BEDRecord("chr2", 500, 10010),
        BEDRecord("chr3", 1, 10),
        BEDRecord("chr3", 100, 200),
    ]


def test_sort_bed_by_bamfile__unknown_contigs_2():
    handle: Any = Mock()
    handle.references = []
    regions = list(_UNSORTED_RECORDS)

    sort_bed_by_bamfile(handle, regions)
    assert regions == [
        BEDRecord("chr1", 0, 1234),
        BEDRecord("chr2", 500, 10010),
        BEDRecord("chr3", 1, 10),
        BEDRecord("chr3", 100, 200),
    ]
