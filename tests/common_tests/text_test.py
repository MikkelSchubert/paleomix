# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

from collections.abc import Iterable
from typing import AnyStr, Callable, NoReturn, TypeVar

import pytest

from paleomix.common import text
from paleomix.common.text import format_timespan

T = TypeVar("T")

###############################################################################
###############################################################################
# Tests for 'padded_table'


def test_format_timespan__seconds() -> None:
    assert format_timespan(0) == "0.0s"
    assert format_timespan(1.2234) == "1.2s"
    assert format_timespan(59.9) == "59.9s"
    assert format_timespan(60) == "1:00s"
    assert format_timespan(27 * 60 + 30.9) == "27:31s"
    assert format_timespan(59 * 60 + 59) == "59:59s"
    assert format_timespan(1 * 60 * 60) == "1:00:00s"
    assert format_timespan(15 * 60 * 60 + 38 * 60 + 17) == "15:38:17s"
    assert format_timespan(123 * 60 * 60 + 38 * 60 + 17) == "123:38:17s"


###############################################################################
###############################################################################
# Tests for 'padded_table'


def padded_table(
    table: Iterable[str | Iterable[object]],
    min_padding: int = 4,
) -> list[str]:
    return list(text.padded_table(table, min_padding=min_padding))


def test_padded_table__empty() -> None:
    assert padded_table(()) == []


def test_padded_table__single_line() -> None:
    table = [(1, 20, 3000)]
    expected = ["1    20    3000"]
    assert expected == padded_table(table)


def test_padded_table__two_lines() -> None:
    table = [(1, 20, 3000), (3000, 20, 1)]
    expected = ["1       20    3000", "3000    20    1"]
    assert expected == padded_table(table)


def test_padded_table__three_lines() -> None:
    table = [(1, 20, 3000), (3000, 20, 1), (1, 2, 30)]
    expected = ["1       20    3000", "3000    20    1", "1       2     30"]
    assert expected == padded_table(table)


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


@pytest.mark.parametrize(("table", "expected"), _PT_PERMUTATIONS)
def test_padded_table__with_text(
    table: tuple[list[object], ...], expected: str
) -> None:
    assert expected == padded_table(table)


###############################################################################
###############################################################################
# Tests for 'parse_linse'


def _this(value: str, length: int) -> tuple[str, int]:
    return (value, length)


def parse_lines(
    lines: Iterable[AnyStr],
    parser: Callable[[AnyStr, int], T],
) -> list[T]:
    return list(text.parse_lines(lines, parser))


def test_parse_lines__empty_file() -> None:
    def _assert_false(_value: str, _length: int) -> NoReturn:
        pytest.fail("parse should not be called")  # pragma: no coverage

    value: list[str] = []
    assert parse_lines(value, _assert_false) == []


def test_parse_lines__single() -> None:
    assert parse_lines(["abc line1 \n"], _this) == [("abc line1", 9)]


_PT_COMMENTS_AND_EMPTY_LINES = (
    ["# comment\n", "abc line1 \n", "def line2 \n"],
    ["abc line1 \n", " # comment\n", "def line2 \n"],
    ["abc line1 \n", "def line2 \n", "   # comment\n"],
    ["\n", "abc line1 \n", "def line2 \n"],
    ["abc line1 \n", " \n", "def line2 \n"],
    ["abc line1 \n", "def line2 \n", "   \n"],
)


@pytest.mark.parametrize("lines", _PT_COMMENTS_AND_EMPTY_LINES)
def test_parse_lines__comments_and_empty_lines(lines: list[str]) -> None:
    expected = [("abc line1", 9), ("def line2", 9)]
    assert parse_lines(lines, _this) == expected


@pytest.mark.parametrize("postfix", ["\r", "\n", "\r\n"])
def test_parse_lines__padding__newlines(postfix: str) -> None:
    expected = [("abc line1", 9), ("def line2", 9)]

    line_1 = "abc line1 " + postfix
    line_2 = "def line2 " + postfix
    assert expected == parse_lines([line_1, line_2], _this)


def test_parse_lines__uncallable() -> None:
    value: list[str] = []
    with pytest.raises(TypeError):
        parse_lines(value, 1)  # pyright: ignore[reportArgumentType]


def test_parse_lines__binary() -> None:
    lines = [b"# foo", b"12", b"3456"]
    expected = [(12, 2), (3456, 4)]

    def parser(value: AnyStr, length: int) -> tuple[int, int]:
        return (int(value), length)

    assert parse_lines(lines, parser) == expected
