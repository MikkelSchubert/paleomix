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

import io
from pathlib import Path

import pytest

from paleomix.tools.validate_fasta import FormatError, check_fasta_file, main


def test_empty_file() -> None:
    with pytest.raises(FormatError, match="File does not contain any sequences"):
        check_fasta_file(io.BytesIO())


def test_simple_sequences() -> None:
    handle = io.BytesIO(b">foo\nACGT\nAA\n>bar\nACtGATA\n>zod\nAcct\nGTAC\nT")

    check_fasta_file(handle)


def test_simple_sequences_carriage_returns() -> None:
    handle = io.BytesIO(b">foo\r\nACGT\r\nAA\r\n>bar\r\nACtGATA\r\n")

    with pytest.raises(FormatError, match="File contains carriage-returns"):
        check_fasta_file(handle)


def test_unnamed_sequence_1() -> None:
    handle = io.BytesIO(b">\nACGT")

    with pytest.raises(FormatError, match="Unnamed fasta sequence"):
        check_fasta_file(handle)


def test_unnamed_sequence_2() -> None:
    handle = io.BytesIO(b"> foo\nACGT")

    with pytest.raises(FormatError, match="Unnamed fasta sequence"):
        check_fasta_file(handle)


def test_invalid_sequence_name() -> None:
    handle = io.BytesIO(b">=bar\nACGT\n")

    with pytest.raises(FormatError, match="Invalid name '=bar'"):
        check_fasta_file(handle)


def test_duplicate_sequence_name() -> None:
    handle = io.BytesIO(b">foo\nACGT\n>foo\nAACC")

    with pytest.raises(FormatError, match="Duplicate sequence name 'foo'"):
        check_fasta_file(handle)


def test_invalid_sequence() -> None:
    handle = io.BytesIO(b">foo\nACxGT")

    with pytest.raises(FormatError, match="Invalid characters 'x' in sequence"):
        check_fasta_file(handle)


def test_empty_sequence_1() -> None:
    handle = io.BytesIO(b">foo\n>bar\nACGT")

    with pytest.raises(FormatError, match="Empty sequence found"):
        check_fasta_file(handle)


def test_empty_sequence_2() -> None:
    handle = io.BytesIO(b">foo\nACGT\n>bar")

    with pytest.raises(FormatError, match="File ends with an empty sequence"):
        check_fasta_file(handle)


def test_no_sequence_header() -> None:
    handle = io.BytesIO(b"ACGT")

    with pytest.raises(FormatError, match="Expected FASTA header, found 'ACGT'"):
        check_fasta_file(handle)


def test_empty_lines_in_sequence_1() -> None:
    handle = io.BytesIO(b">foo\nACGT\n\nAA\n")

    with pytest.raises(FormatError, match="Empty lines not allowed in sequences"):
        check_fasta_file(handle)


def test_empty_lines_in_sequence_2() -> None:
    handle = io.BytesIO(b">foo\n\nGTTG\nAA\n")

    with pytest.raises(FormatError, match="Expected FASTA sequence, found empty line"):
        check_fasta_file(handle)


def test_whitespace() -> None:
    handle = io.BytesIO(b">foo\nACGT\nAA\n\n\n>bar\nACtGATA\n\n")

    check_fasta_file(handle)


def test_main_ok(tmp_path: Path) -> None:
    filepath = tmp_path / "file.fasta"
    filepath.write_text(">foo\nACGTA")

    assert main([str(filepath)]) == 0


def test_main_err(tmp_path: Path) -> None:
    filepath = tmp_path / "file.fasta"
    filepath.write_text(">foo\nA\nCGTA")

    assert main([str(filepath)]) == 1
