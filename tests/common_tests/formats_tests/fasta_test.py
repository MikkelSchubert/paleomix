import StringIO
import random
import nose.tools

from pypeline.common.formats.fasta import *

def assert_list_equals(aa, bb):
    aa = list(aa)
    bb = list(bb)

    assert aa == bb



_SEQ_FRAG = "AAGTCC" # len() = 6


################################################################################
################################################################################
## Tests for wrap_fasta

def test_wrap_fasta__partial_line_test():
    expected = ">foobar\n%s\n" % (_SEQ_FRAG, )
    result = wrap_fasta("foobar", _SEQ_FRAG)
    assert result == expected
    
def test_wrap_fasta__complete_line_test():
    expected = ">barfoo\n%s\n" % (_SEQ_FRAG * 10, )
    result = wrap_fasta("barfoo", _SEQ_FRAG * 10)
    assert result == expected

def test_wrap_fasta__multiple_lines():
    expected = ">foobar\n%s\n%s\n" \
        % (_SEQ_FRAG * 10, _SEQ_FRAG * 5)
    result = wrap_fasta("foobar", _SEQ_FRAG * 15)
    assert result == expected




################################################################################
################################################################################
## Tests for print_fasta

def test_print_fasta__multiple_lines():
    expected = ">foobar\n%s\n%s\n" \
        % (_SEQ_FRAG * 10, _SEQ_FRAG * 5)
    stringf = StringIO.StringIO()
    print_fasta("foobar", _SEQ_FRAG * 15, stringf)
    assert stringf.getvalue() == expected




################################################################################
################################################################################
## Tests for parse_fasta


def _get_sequence(length):
    return "".join(random.choice("ACGT") for i in range(length))

def _get_lines( name, sequence):
    lines = [">" + name]
    while sequence:
        lines.append(sequence[:60] + "\n")
        sequence = sequence[60:]
    return lines

        

def test_parse_fasta__no_records():
    assert_list_equals(parse_fasta([]), [])


def test_parse_fasta__single_record():
    def test_function(length):
        name     = "foobar%i" % length
        sequence = _get_sequence(length)
        lines    = _get_lines(name, sequence)
        
        assert_list_equals(parse_fasta(lines), [(name, sequence)])

    for length in (45, 60, 90, 120):
        yield test_function, length


def test_parse_fasta__n_records():
    def test_function(nrecords):
        lines, records = [], []
        for _ in range(nrecords):
            length   = random.randint(5, 240)
            sequence = _get_sequence(length)
            name     = "random%i" % length
                
            records.append((name, sequence))
            lines.extend(_get_lines(name, sequence))
    
        assert_list_equals(parse_fasta(lines), records)

    for nrecords in range(2, 5):
        yield test_function, nrecords


@nose.tools.raises(ValueError)
def test_parse_fasta__empty_record_name_only__nothing_else():
    list(parse_fasta([">fasta1\n"]))

@nose.tools.raises(ValueError)
def test_parse_fasta__empty_record_name_only__first():
    list(parse_fasta([">fasta1\n", ">fasta2\n", "AGTC\n"]))

@nose.tools.raises(ValueError)
def test_parse_fasta__empty_record__middle():
    lines = [">fasta0\n", "ACGT\n", ">fasta1\n", ">fasta2\n", "AGTC\n"]
    list(parse_fasta(lines))

@nose.tools.raises(ValueError)
def test_parse_empty_record_last():
    lines = [">fasta1\n", "ACGT\n", ">fasta2\n"]
    list(parse_fasta(lines))
    
@nose.tools.raises(ValueError)
def test_parse_fasta__missing_name__alone():
    lines = ["ACGT\n"]
    list(parse_fasta(lines))

@nose.tools.raises(ValueError)
def test_parse_fasta__missing_name__with_others():
    lines = ["ACGT\n", ">Foo\n", "ACGGTA\n"]
    list(parse_fasta(lines))
    



################################################################################
################################################################################
## Tests for 'read_fasta'

def test_read_fasta__uncompressed():
    expected = [("This is FASTA!", "ACGTN"),
                ("This is ALSO FASTA!", "CGTNA")]
    results  = list(read_fasta("tests/data/fasta_file.fasta"))

    assert results == expected

def test_read_fasta__compressed():
    expected = [("This is GZipped FASTA!", "ACGTN"),
                ("This is ALSO GZipped FASTA!", "CGTNA")]
    results  = list(read_fasta("tests/data/fasta_file.fasta.gz"))

    assert results == expected
