#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import itertools
from typing import Any, AnyStr, Callable, Iterable, Iterator, TypeVar

T = TypeVar("T")


def format_timespan(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds // 60:.0f}:{seconds % 60:02.0f}s"
    else:
        return "{:.0f}:{:02.0f}:{:02.0f}s".format(
            seconds // 3600,
            (seconds % 3600) // 60,
            seconds % 60,
        )


def padded_table(
    table: Iterable[str | Iterable[Any]],
    min_padding: int = 4,
) -> Iterator[str]:
    """Takes a sequence of iterables, each of which represents a row in a
    table. Values are converted to string, and padded with whitespace such that
    each column is separated from its adjacent columns by at least 4 spaces.
    Empty cells or whitespace in values are not allowed.

    If a string is included instead of a row, this value is added as is. Note
    that these lines should be whitespace only, or start with a '#' if the
    resulting table is to be readable with 'parse_padded_table'.
    """
    str_rows: list[str | list[str]] = []
    max_sizes: list[int] = []
    for row in table:
        if not isinstance(row, str):
            row = list(map(str, row))
            row_sizes = list(map(len, row))
            max_sizes = list(
                map(max, itertools.zip_longest(max_sizes, row_sizes, fillvalue=0))
            )

        str_rows.append(row)

    sizes = [(size + min_padding) for size in max_sizes]
    for row in str_rows:
        if not isinstance(row, str):
            row = "".join(
                field.ljust(padding) for (field, padding) in zip(row, sizes)
            ).rstrip()
        yield row


def parse_lines(
    lines: Iterable[AnyStr],
    parser: Callable[[AnyStr, int], T],
) -> Iterator[T]:
    """Parses a set of lines using the supplied callable:
        lambda (line, length): ...

    Supports the parser functions available in 'pysam': asGTF, asBED, etc.
    """
    if not callable(parser):
        raise TypeError(f"'parser' must be a callable, not {parser!r}")

    for line in lines:
        stripped = line.lstrip()
        if stripped and stripped[0] not in ("#", 35):
            stripped = line.rstrip()
            yield parser(stripped, len(stripped))
