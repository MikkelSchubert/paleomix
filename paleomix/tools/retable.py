#!/usr/bin/python
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
"""
Reads and pretty-prints whitespace separated tabular data.
"""
import sys
import argparse
import fileinput

from paleomix.common.text import padded_table


def split_line(line):
    sections = line.rstrip('\r\n').split('#', 1)
    if len(sections) == 1:
        sections.append(None)

    text, comment = sections
    columns = text.split()

    if comment is not None:
        if not columns:
            return '#' + comment

        columns.append('#' + comment)

    return columns


def parse_args(argv):
    parser = argparse.ArgumentParser(prog="paleomix retable")
    parser.add_argument("file", nargs="*",
                        help="One or more text files containing tabular data.")

    return parser.parse_args(argv)


def main(argv):
    """Main function; takes a list of arguments but excluding sys.argv[0]."""
    args = parse_args(argv)

    rows = []
    for line in fileinput.input(args.file):
        rows.append(split_line(line))

    for line in padded_table(rows):
        print(line)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
