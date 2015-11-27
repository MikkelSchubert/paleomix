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
# Disable warnings on strange function names
# pylint: disable=C0103


import collections

import nose.tools
from nose.tools import assert_equal

from paleomix.common.text import \
    TableError, \
    padded_table, \
    parse_padded_table, \
    parse_lines, \
    parse_lines_by_contig


###############################################################################
###############################################################################
# Tests for 'padded_table'

def _padded_table(*args, **kwargs):
    return list(padded_table(*args, **kwargs))


def test_padded_table__empty():
    assert_equal(_padded_table(()), [])


def test_padded_table__single_line():
    table = [(1, 20, 3000)]
    expected = ["1    20    3000"]
    assert_equal(expected, _padded_table(table))


def test_padded_table__two_lines():
    table = [(1, 20, 3000),
             (3000, 20, 1)]
    expected = ["1       20    3000",
                "3000    20    1"]
    assert_equal(expected, _padded_table(table))


def test_padded_table__three_lines():
    table = [(1, 20, 3000),
             (3000, 20, 1),
             (1, 2, 30)]
    expected = ["1       20    3000",
                "3000    20    1",
                "1       2     30"]
    assert_equal(expected, _padded_table(table))


def test_padded_table__with_text():
    row_1, line_1 = (1, 20, 3000), "1       20    3000"
    row_2, line_2 = (3000, 20, 1), "3000    20    1"
    comment = "# An insightful comment goes here"

    def _do_test_padded_table__padding__with_text(table, expected):
        assert_equal(expected, _padded_table(table))
    yield _do_test_padded_table__padding__with_text, \
        [comment, row_1, row_2], [comment, line_1, line_2]
    yield _do_test_padded_table__padding__with_text, \
        [row_1, comment, row_2], [line_1, comment, line_2]
    yield _do_test_padded_table__padding__with_text, \
        [row_1, row_2, comment], [line_1, line_2, comment]


@nose.tools.raises(TableError)
def test_padded_table__misshapen_table_1():
    _padded_table([(1, 20, 3000),
                   (3000, 20),
                   (1, 2, 3)])


@nose.tools.raises(TableError)
def test_padded_table__misshapen_table_2():
    _padded_table([(1, 20, 3000),
                   (3000, 20, 1, 0),
                   (1, 2, 3)])


###############################################################################
###############################################################################
# Tests for 'parse_padded_table'

def _parse_padded_table(*args, **kwargs):
    return list(parse_padded_table(*args, **kwargs))


def test_parse_padded_table__empty():
    assert_equal([], _parse_padded_table([]))


def test_parse_padded_table__header_only():
    assert_equal([], _parse_padded_table(["A  B  C  D"]))


def test_parse_padded_table__single_row():
    table = ["A    B    C    D",
             "4    3    2    1"]
    expected = [{"A": '4', "B": '3', "C": '2', "D": '1'}]
    assert_equal(expected, _parse_padded_table(table))


def test_parse_padded_table__two_rows():
    table = ["A     B     C     D",
             "4     3     2     1",
             "AB    CD    EF    GH"]
    expected = [{"A": '4', "B": '3', "C": '2', "D": '1'},
                {"A": "AB", "B": "CD", "C": "EF", "D": "GH"}]
    assert_equal(expected, _parse_padded_table(table))


# Any amount of whitespace is allowed
def test_parse_padded_table__single_row__with_whitespace():
    table = ["A   B    C       E F",
             "1        0  1    2   3"]
    expected = [{"A": '1', "B": '0', "C": '1', "E": '2', "F": '3'}]
    assert_equal(expected, _parse_padded_table(table))


# Other whitespace should be ignored
def test_parse_padded_table__single_row__with_tabs():
    table = ["A\t\t\t\tB",
             "1\t\t\t\t0"]
    expected = [{"A": '1', "B": '0'}]
    assert_equal(expected, _parse_padded_table(table))


def test_padded_table__comments_and_empty_lines():
    def _do_test_padded_table__comments(lines):
        expected = [{"A": '3000', "B": '20', "C": '1'}]
        assert_equal(expected, _parse_padded_table(lines))
    line_1 = "A         B      C"
    line_2 = "3000      20      1"

    yield _do_test_padded_table__comments, ["# comment", line_1, line_2]
    yield _do_test_padded_table__comments, [line_1, " # comment", line_2]
    yield _do_test_padded_table__comments, [line_1, line_2, "  # comment"]

    yield _do_test_padded_table__comments, ["", line_1, line_2]
    yield _do_test_padded_table__comments, [line_1, "  ", line_2]
    yield _do_test_padded_table__comments, [line_1, line_2, "   "]


def test_padded_table__newlines():
    expected = [{"A": '3000', "B": '20', "C": '1'}]

    def _do_test_padded_table__padding__comments(postfix):
        line_1 = "A         B       C" + postfix
        line_2 = "3000      20      1" + postfix
        assert_equal(expected, _parse_padded_table([line_1, line_2]))

    yield _do_test_padded_table__padding__comments, "\r"
    yield _do_test_padded_table__padding__comments, "\n"
    yield _do_test_padded_table__padding__comments, "\r\n"


def test_padded_table__padding__comments__whitespace():
    expected = [{"A": '3000', "B": '20', "C": '1'}]
    lines = ["A         B       C",
             "3000      20      1",
             "  # useless comment"]
    assert_equal(expected, _parse_padded_table(lines))


@nose.tools.raises(TableError)
def test_parse_padded_table__malformed_table_0():
    table = ["A    B    C    D",
             "4    3    2"]
    _parse_padded_table(table)


@nose.tools.raises(TableError)
def test_parse_padded_table__malformed_table_1():
    table = ["A    B    C    D",
             "4    3    2    1    0"]
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
    assert_equal(_parse_lines([], _assert_false), [])


def test_parse_lines__single():
    assert_equal(_parse_lines(["abc line1 \n"], _this), [("abc line1", 9)])


def test_parse_lines__comments_and_empty_lines():
    def _do_test_parse_lines__comments(lines):
        expected = [("abc line1", 9), ("def line2", 9)]
        assert_equal(_parse_lines(lines, _this), expected)
    yield _do_test_parse_lines__comments, \
        ["# comment\n", "abc line1 \n", "def line2 \n"]
    yield _do_test_parse_lines__comments, \
        ["abc line1 \n", " # comment\n", "def line2 \n"]
    yield _do_test_parse_lines__comments, \
        ["abc line1 \n", "def line2 \n", "   # comment\n"]

    yield _do_test_parse_lines__comments, \
        ["\n", "abc line1 \n", "def line2 \n"]
    yield _do_test_parse_lines__comments, \
        ["abc line1 \n", " \n", "def line2 \n"]
    yield _do_test_parse_lines__comments, \
        ["abc line1 \n", "def line2 \n", "   \n"]


def test_parse_lines__padding__newlines():
    expected = [("abc line1", 9), ("def line2", 9)]

    def _do_test_padded_table__padding__comments(postfix):
        line_1 = "abc line1 " + postfix
        line_2 = "def line2 " + postfix
        assert_equal(expected, _parse_lines([line_1, line_2], _this))

    yield _do_test_padded_table__padding__comments, "\r"
    yield _do_test_padded_table__padding__comments, "\n"
    yield _do_test_padded_table__padding__comments, "\r\n"


@nose.tools.raises(TypeError)
def test_parse_lines__uncallable():
    _parse_lines([], 1)


###############################################################################
###############################################################################
# Tests for 'parse_lines_by_contig'
_RecordMock = collections.namedtuple("_RecordMock", "contig value")


def test_parse_lines_by_contig__single_contig():
    lines = ["abc line1 \n", "abc line2 \n"]

    def _parse(line, length):
        assert_equal(len(line), length)
        return _RecordMock(*line.split())

    expected = {"abc": [_RecordMock("abc", "line1"),
                        _RecordMock("abc", "line2")]}
    assert_equal(parse_lines_by_contig(lines, _parse), expected)


def test_parse_lines__two_contigs():
    lines = ["abc line1 \n", "def line2 \n"]

    def _parse(line, length):
        assert_equal(len(line), length)
        return _RecordMock(*line.split())

    expected = {"abc": [_RecordMock("abc", "line1")],
                "def": [_RecordMock("def", "line2")]}
    assert_equal(parse_lines_by_contig(lines, _parse), expected)
