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
import re
import types
import collections


class TableError(RuntimeError):
    pass


_MIN_PADDING = 4
_WHITESPACE_OR_EMPTY = re.compile(r"\s|^$")
def padded_table(table):
    """Takes a sequence of iterables, each of which represents a row in a table.
    Values are converted to string, and padded with whitespace such that each
    column is separated from its adjacent columns by at least 4 spaces. Empty
    cells or whitespace in values are not allowed.

    If a string is included instead of a row, this value is added as is. Note
    that these lines should be whitespace only, or start with a '#' if the
    resulting table is to be readable with 'parse_padded_table'."""
    str_rows  = []
    nsizes, sizes = None, []
    for row in table:
        if not isinstance(row, types.StringTypes):
            row = map(str, row)
            if (len(row) != nsizes):
                if nsizes is not None:
                    raise TableError("Malformed table; rows with different number of columns: %r" % row)
                nsizes = len(row)
                sizes  = [0] * nsizes
            sizes = map(max, zip(sizes, map(len, row)))
        str_rows.append(row)

    sizes = [(size + _MIN_PADDING) for size in sizes]
    for row in str_rows:
        if not isinstance(row, types.StringTypes):
            row = "".join(field.ljust(padding) for (field, padding) in zip(row, sizes)).rstrip() # pragma: no coverage
        yield row


_SPLIT_BY_WHITESPACE = re.compile((" " * _MIN_PADDING) + "+")
def parse_padded_table(lines, header = None):
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        elif header is None:
            header  = _SPLIT_BY_WHITESPACE.split(stripped)
            nheader = len(header)
            continue

        fields = _SPLIT_BY_WHITESPACE.split(stripped)
        if len(fields) != nheader:
            raise TableError("Malformed table; #columns does not match header: %r vs %r" % (header, fields))

        yield dict(zip(header, fields))


def parse_lines(lines, parser):
    """Parses a set of lines using the supplied callable: lambda (line, length): ...
    Supports the parser functions available in 'pysam': asGTF, asBED, etc."""
    if not isinstance(parser, collections.Callable):
        raise TypeError("'parser' must be a callable, not %r" % parser.__class__.__name__)

    for line in lines:
        stripped = line.lstrip()
        if stripped and not stripped.startswith("#"):
            stripped = line.rstrip()
            yield parser(stripped, len(stripped))


def parse_lines_by_contig(lines, parser):
    """Reads the lines of a text file, parsing each line with the specified
    parser, and aggregating results by the 'contig' property of reach record."""
    table = {}
    for record in parse_lines(lines, parser):
        try:
            table[record.contig].append(record)
        except KeyError:
            table[record.contig] = [record]

    return table
