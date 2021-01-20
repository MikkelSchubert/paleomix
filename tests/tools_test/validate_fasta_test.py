import io

import pytest

from paleomix.tools.validate_fasta import main, check_fasta_file, FormatError


def test_empty_file():
    with pytest.raises(FormatError, match="File does not contain any sequences"):
        check_fasta_file(io.BytesIO())


def test_simple_sequences():
    handle = io.BytesIO(b">foo\nACGT\nAA\n>bar\nACtGATA\n>zod\nAcct\nGTAC\nT")

    check_fasta_file(handle)


def test_simple_sequences_carriage_returns():
    handle = io.BytesIO(b">foo\r\nACGT\r\nAA\r\n>bar\r\nACtGATA\r\n")

    with pytest.raises(FormatError, match="File contains carriage-returns"):
        check_fasta_file(handle)


def test_unnamed_sequence_1():
    handle = io.BytesIO(b">\nACGT")

    with pytest.raises(FormatError, match="Unnamed fasta sequence"):
        check_fasta_file(handle)


def test_unnamed_sequence_2():
    handle = io.BytesIO(b"> foo\nACGT")

    with pytest.raises(FormatError, match="Unnamed fasta sequence"):
        check_fasta_file(handle)


def test_invalid_sequence_name():
    handle = io.BytesIO(b">=bar\nACGT\n")

    with pytest.raises(FormatError, match="Invalid name '=bar'"):
        check_fasta_file(handle)


def test_duplicate_sequence_name():
    handle = io.BytesIO(b">foo\nACGT\n>foo\nAACC")

    with pytest.raises(FormatError, match="Duplicate sequence name 'foo'"):
        check_fasta_file(handle)


def test_invalid_sequence():
    handle = io.BytesIO(b">foo\nACxGT")

    with pytest.raises(FormatError, match="Invalid characters 'x' in sequence"):
        check_fasta_file(handle)


def test_empty_sequence_1():
    handle = io.BytesIO(b">foo\n>bar\nACGT")

    with pytest.raises(FormatError, match="Empty sequence found"):
        check_fasta_file(handle)


def test_empty_sequence_2():
    handle = io.BytesIO(b">foo\nACGT\n>bar")

    with pytest.raises(FormatError, match="File ends with an empty sequence"):
        check_fasta_file(handle)


def test_no_sequence_header():
    handle = io.BytesIO(b"ACGT")

    with pytest.raises(FormatError, match="Expected FASTA header, found 'ACGT'"):
        check_fasta_file(handle)


def test_empty_lines_in_sequence_1():
    handle = io.BytesIO(b">foo\nACGT\n\nAA\n")

    with pytest.raises(FormatError, match="Empty lines not allowed in sequences"):
        check_fasta_file(handle)


def test_empty_lines_in_sequence_2():
    handle = io.BytesIO(b">foo\n\nGTTG\nAA\n")

    with pytest.raises(FormatError, match="Expected FASTA sequence, found empty line"):
        check_fasta_file(handle)


def test_whitespace():
    handle = io.BytesIO(b">foo\nACGT\nAA\n\n\n>bar\nACtGATA\n\n")

    check_fasta_file(handle)


def test_main_ok(tmp_path):
    filepath = tmp_path / "file.fasta"
    filepath.write_text(">foo\nACGTA")

    assert main(str(filepath)) == 0


def test_main_err(tmp_path):
    filepath = tmp_path / "file.fasta"
    filepath.write_text(">foo\nA\nCGTA")

    assert main(str(filepath)) == 1
