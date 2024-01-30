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

import re
import sys
from pathlib import Path
from typing import IO, TYPE_CHECKING, NoReturn

from paleomix.common.argparse import ArgumentParser
from paleomix.common.formats import FormatError

if TYPE_CHECKING:
    from paleomix.common import argparse

# Standard nucleotides + UIPAC codes
_VALID_CHARS_STR = b"ACGTNRYSWKMBDHV"
_VALID_CHARS = frozenset(_VALID_CHARS_STR.upper() + _VALID_CHARS_STR.lower())
_NA, _IN_HEADER, _IN_SEQUENCE, _IN_WHITESPACE = range(4)

_SEQUENCE_NAME = re.compile(b"[!-()+-<>-~][!-~]*")


def _raise_error(line: str, value: bytes | None = None) -> NoReturn:
    if value is not None:
        line = line.format(value.decode("utf-8", errors="replace"))

    raise FormatError(line)


def _validate_fasta_header(line: bytes, sequence_names: set[bytes]) -> None:
    name = line.split(b" ", 1)[0][1:]
    if not name:
        _raise_error("Unnamed fasta sequence")
    elif not _SEQUENCE_NAME.match(name):
        _raise_error("Invalid name {!r}", name)
    elif name in sequence_names:
        _raise_error("Duplicate sequence name {!r}", name)

    sequence_names.add(name)


def _validate_fasta_line(
    line: bytes, valid_chars: frozenset[int] = _VALID_CHARS
) -> None:
    invalid_chars = bytes(frozenset(line) - valid_chars)

    if invalid_chars:
        _raise_error("Invalid characters {!r} in sequence", invalid_chars)


def check_fasta_file(handle: IO[bytes]) -> None:
    linenum = 0
    sequence_names: set[bytes] = set()
    state, linelength, linelengthchanged = _NA, -1, False
    try:
        for linenum, line in enumerate(handle, start=1):
            # Only \n is allowed as not all tools handle \r
            line = line.rstrip(b"\n")

            if line.endswith(b"\r"):
                _raise_error(
                    "File contains carriage-returns ('\\r')! Please convert file to "
                    "unix format, using e.g. dos2unix."
                )
            elif line.startswith(b">"):
                if state in (_NA, _IN_SEQUENCE, _IN_WHITESPACE):
                    _validate_fasta_header(line, sequence_names)
                    state = _IN_HEADER
                    linelength = -1
                    linelengthchanged = False
                else:
                    assert state == _IN_HEADER, state
                    _raise_error("Empty sequence found")
            elif line:
                if state == _IN_SEQUENCE:
                    _validate_fasta_line(line)
                    # Only the final line in a sequence is allowed to differ in length
                    # from the proceeding lines. Otherwise the file cannot be faidx'd.
                    if linelengthchanged or linelength < len(line):
                        _raise_error("Lines in FASTQ files must be of same length")
                    elif linelength != len(line):
                        linelengthchanged = True
                elif state == _IN_HEADER:
                    _validate_fasta_line(line)
                    linelength = len(line)
                    state = _IN_SEQUENCE
                elif state == _NA:
                    _raise_error("Expected FASTA header, found {!r}", line)
                else:
                    assert state == _IN_WHITESPACE, state
                    _raise_error("Empty lines not allowed in sequences")
            elif state in (_NA, _IN_WHITESPACE):
                continue
            elif state == _IN_HEADER:
                _raise_error("Expected FASTA sequence, found empty line")
            else:
                assert state == _IN_SEQUENCE, state
                state = _IN_WHITESPACE
    except FormatError as error:
        _raise_error(f"Line {linenum}: {error}")

    if state == _NA:
        _raise_error("File does not contain any sequences")
    elif state == _IN_HEADER:
        _raise_error("File ends with an empty sequence")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = ArgumentParser("paleomix :validate_fasta")
    parser.add_argument("fasta", type=Path)

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    with args.fasta.open("rb") as handle:
        try:
            check_fasta_file(handle)
        except FormatError as error:
            print(error)
            return 1

    print("No errors found.")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv[1:]))
