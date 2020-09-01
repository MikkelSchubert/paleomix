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
import collections

import pytest

from paleomix.common.text import (
    TableError,
    padded_table,
    parse_lines,
    parse_lines_by_contig,
    parse_padded_table,
)

###############################################################################
###############################################################################
# Tests for 'padded_table'


def _padded_table(*args, **kwargs):
    return list(padded_table(*args, **kwargs))


def test_padded_table__empty():
    assert _padded_table(()) == []


def test_padded_table__single_line():
    table = [(1, 20, 3000)]
    expected = ["1    20    3000"]
    assert expected == _padded_table(table)


def test_padded_table__two_lines():
    table = [(1, 20, 3000), (3000, 20, 1)]
    expected = ["1       20    3000", "3000    20    1"]
    assert expected == _padded_table(table)


def test_padded_table__three_lines():
    table = [(1, 20, 3000), (3000, 20, 1), (1, 2, 30)]
    expected = ["1       20    3000", "3000    20    1", "1       2     30"]
    assert expected == _padded_table(table)


_PT_COMMENT = "# An insightful comment goes here"
_PT_ROW_1 = (1, 20, 3000)
_PT_ROW_2 = (3000, 20, 1)
_PT_LINE_1 = "1       20    3000"
_PT_LINE_2 = "3000    20    1"

_PT_PERMUTATIONS = (
    ([_PT_COMMENT, _PT_ROW_1, _PT_ROW_2], [_PT_COMMENT, _PT_LINE_1, _PT_LINE_2]),
    ([_PT_ROW_1, _PT_COMMENT, _PT_ROW_2], [_PT_LINE_1, _PT_COMMENT, _PT_LINE_2]),
    ([_PT_ROW_1, _PT_ROW_2, _PT_COMMENT], [_PT_LINE_1, _PT_LINE_2, _PT_COMMENT]),
)


@pytest.mark.parametrize("table, expected", _PT_PERMUTATIONS)
def test_padded_table__with_text(table, expected):
    assert expected == _padded_table(table)


###############################################################################
###############################################################################
# Tests for 'parse_padded_table'


def _parse_padded_table(*args, **kwargs):
    return list(parse_padded_table(*args, **kwargs))


def test_parse_padded_table__empty():
    assert [] == _parse_padded_table([])


def test_parse_padded_table__header_only():
    assert [] == _parse_padded_table(["A  B  C  D"])


def test_parse_padded_table__single_row():
    table = ["A    B    C    D", "4    3    2    1"]
    expected = [{"A": "4", "B": "3", "C": "2", "D": "1"}]
    assert expected == _parse_padded_table(table)


def test_parse_padded_table__two_rows():
    table = ["A     B     C     D", "4     3     2     1", "AB    CD    EF    GH"]
    expected = [
        {"A": "4", "B": "3", "C": "2", "D": "1"},
        {"A": "AB", "B": "CD", "C": "EF", "D": "GH"},
    ]
    assert expected == _parse_padded_table(table)


# Any amount of whitespace is allowed
def test_parse_padded_table__single_row__with_whitespace():
    table = ["A   B    C       E F", "1        0  1    2   3"]
    expected = [{"A": "1", "B": "0", "C": "1", "E": "2", "F": "3"}]
    assert expected == _parse_padded_table(table)


# Other whitespace should be ignored
def test_parse_padded_table__single_row__with_tabs():
    table = ["A\t\t\t\tB", "1\t\t\t\t0"]
    expected = [{"A": "1", "B": "0"}]
    assert expected == _parse_padded_table(table)


_PT_COMMENTS_LINE_AND_EMPTY_1 = "A         B      C"
_PT_COMMENTS_LINE_AND_EMPTY_2 = "3000      20      1"

__PT_COMMENTS_LINE_AND_EMPTY = (
    ["# comment", _PT_COMMENTS_LINE_AND_EMPTY_1, _PT_COMMENTS_LINE_AND_EMPTY_2],
    [_PT_COMMENTS_LINE_AND_EMPTY_1, " # comment", _PT_COMMENTS_LINE_AND_EMPTY_2],
    [_PT_COMMENTS_LINE_AND_EMPTY_1, _PT_COMMENTS_LINE_AND_EMPTY_2, "  # comment"],
    ["", _PT_COMMENTS_LINE_AND_EMPTY_1, _PT_COMMENTS_LINE_AND_EMPTY_2],
    [_PT_COMMENTS_LINE_AND_EMPTY_1, "  ", _PT_COMMENTS_LINE_AND_EMPTY_2],
    [_PT_COMMENTS_LINE_AND_EMPTY_1, _PT_COMMENTS_LINE_AND_EMPTY_2, "   "],
)


@pytest.mark.parametrize("lines", __PT_COMMENTS_LINE_AND_EMPTY)
def test_padded_table__comments_and_empty_lines(lines):
    expected = [{"A": "3000", "B": "20", "C": "1"}]
    assert expected == _parse_padded_table(lines)


@pytest.mark.parametrize("postfix", ("\r", "\n", "\r\n"))
def test_padded_table__newlines(postfix):
    expected = [{"A": "3000", "B": "20", "C": "1"}]

    line_1 = "A         B       C" + postfix
    line_2 = "3000      20      1" + postfix
    assert expected == _parse_padded_table([line_1, line_2])


def test_padded_table__padding__comments__whitespace():
    expected = [{"A": "3000", "B": "20", "C": "1"}]
    lines = ["A         B       C", "3000      20      1", "  # useless comment"]
    assert expected == _parse_padded_table(lines)


def test_parse_padded_table__malformed_table_0():
    table = ["A    B    C    D", "4    3    2"]
    with pytest.raises(TableError):
        _parse_padded_table(table)


def test_parse_padded_table__malformed_table_1():
    table = ["A    B    C    D", "4    3    2    1    0"]
    with pytest.raises(TableError):
        _parse_padded_table(table)


###############################################################################
###############################################################################
# Tests for 'parse_linse'
def _this(*args):
    return args


def _parse_lines(*args, **kwargs):
    return list(parse_lines(*args, **kwargs))


def test_parse_lines__empty_file():
    def _assert_false():
        assert False  # pragma: no coverage

    assert _parse_lines([], _assert_false) == []


def test_parse_lines__single():
    assert _parse_lines(["abc line1 \n"], _this) == [("abc line1", 9)]


_PT_COMMENTS_AND_EMPTY_LINES = (
    ["# comment\n", "abc line1 \n", "def line2 \n"],
    ["abc line1 \n", " # comment\n", "def line2 \n"],
    ["abc line1 \n", "def line2 \n", "   # comment\n"],
    ["\n", "abc line1 \n", "def line2 \n"],
    ["abc line1 \n", " \n", "def line2 \n"],
    ["abc line1 \n", "def line2 \n", "   \n"],
)


@pytest.mark.parametrize("lines", _PT_COMMENTS_AND_EMPTY_LINES)
def test_parse_lines__comments_and_empty_lines(lines):
    expected = [("abc line1", 9), ("def line2", 9)]
    assert _parse_lines(lines, _this) == expected


@pytest.mark.parametrize("postfix", ("\r", "\n", "\r\n"))
def test_parse_lines__padding__newlines(postfix):
    expected = [("abc line1", 9), ("def line2", 9)]

    line_1 = "abc line1 " + postfix
    line_2 = "def line2 " + postfix
    assert expected == _parse_lines([line_1, line_2], _this)


def test_parse_lines__uncallable():
    with pytest.raises(TypeError):
        _parse_lines([], 1)


def parse_lines__binary():
    lines = [b"# foo", b"12", b"3456"]
    expected = [(12, 2), (3456, 4)]

    def parser(value, length):
        return (int(value), length)

    assert parse_lines(lines, parser) == expected


###############################################################################
###############################################################################
# Tests for 'parse_lines_by_contig'
_RecordMock = collections.namedtuple("_RecordMock", "contig value")


def test_parse_lines_by_contig__single_contig():
    lines = ["abc line1 \n", "abc line2 \n"]

    def _parse(line, length):
        assert len(line) == length
        return _RecordMock(*line.split())

    expected = {"abc": [_RecordMock("abc", "line1"), _RecordMock("abc", "line2")]}
    assert parse_lines_by_contig(lines, _parse) == expected


def test_parse_lines__two_contigs():
    lines = ["abc line1 \n", "def line2 \n"]

    def _parse(line, length):
        assert len(line) == length
        return _RecordMock(*line.split())

    expected = {
        "abc": [_RecordMock("abc", "line1")],
        "def": [_RecordMock("def", "line2")],
    }
    assert parse_lines_by_contig(lines, _parse) == expected
