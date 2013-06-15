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
import StringIO
import nose.tools
from nose.tools import assert_equals # pylint: disable=E0611
from pypeline.common.formats.fasta import \
     wrap_fasta, \
     print_fasta, \
     parse_fasta, \
     read_fasta, \
     FASTAError



def assert_list_equals(iter_a, iter_b):
    """Compare two values, after first converting them to lists.
    This enures that lazily generated results can be compared."""
    list_a = list(iter_a)
    list_b = list(iter_b)

    assert_equals(list_a, list_b)



_SEQ_FRAG = "AAGTCC" # len() = 6


################################################################################
################################################################################
## Tests for wrap_fasta

def test_wrap_fasta__partial_line_test():
    expected = ">foobar\n%s\n" % (_SEQ_FRAG, )
    result = wrap_fasta("foobar", _SEQ_FRAG)
    assert_equals(result, expected)

def test_wrap_fasta__complete_line_test():
    expected = ">barfoo\n%s\n" % (_SEQ_FRAG * 10, )
    result = wrap_fasta("barfoo", _SEQ_FRAG * 10)
    assert_equals(result, expected)

def test_wrap_fasta__multiple_lines():
    expected = ">foobar\n%s\n%s\n" \
        % (_SEQ_FRAG * 10, _SEQ_FRAG * 5)
    result = wrap_fasta("foobar", _SEQ_FRAG * 15)
    assert_equals(result, expected)




################################################################################
################################################################################
## Tests for print_fasta


def test_print_fasta__partial_line():
    expected = ">foobar\n%s\n" % (_SEQ_FRAG, )
    stringf = StringIO.StringIO()
    print_fasta("foobar", _SEQ_FRAG, stringf)
    assert_equals(stringf.getvalue(), expected)

def test_print_fasta__complete_line_test():
    expected = ">barfoo\n%s\n" % (_SEQ_FRAG * 10, )
    stringf = StringIO.StringIO()
    print_fasta("barfoo", _SEQ_FRAG * 10, stringf)
    assert_equals(stringf.getvalue(), expected)

def test_print_fasta__multiple_lines():
    expected = ">foobar\n%s\n%s\n" \
        % (_SEQ_FRAG * 10, _SEQ_FRAG * 5)
    stringf = StringIO.StringIO()
    print_fasta("foobar", _SEQ_FRAG * 15, stringf)
    assert_equals(stringf.getvalue(), expected)




################################################################################
################################################################################
## Tests for parse_fasta

def test_parse_fasta__no_records():
    assert_list_equals(parse_fasta([]), [])

def test_parse_fasta__single_record():
    lines    = [">single\n", "TGTTCTCCACCGTGCACAAC\n", "CCTTCATCCA\n"]
    expected = [(("single", None), "TGTTCTCCACCGTGCACAACCCTTCATCCA")]
    assert_list_equals(parse_fasta(lines), expected)

def test_parse_fasta__multiple_records():
    lines    = [">first\n",  "TGTTCTCCACCGTGCACAAC\n", "CCTTCATCCA\n",
                ">Second XT:1:0\n", "GAGAGCTCAGCTAAC\n",
                ">Third\n",  "CGCTGACCAAAAACGGACAG\n", "GGCATTCGGC\n"]
    expected = [(("first", None), "TGTTCTCCACCGTGCACAACCCTTCATCCA"),
                (("Second", "XT:1:0"), "GAGAGCTCAGCTAAC"),
                (("Third", None), "CGCTGACCAAAAACGGACAGGGCATTCGGC")]
    assert_list_equals(parse_fasta(lines), expected)

@nose.tools.raises(FASTAError)
def test_parse_fasta__empty_record_name_only__nothing_else():
    list(parse_fasta([">fasta1\n"]))

@nose.tools.raises(FASTAError)
def test_parse_fasta__empty_record_name_only__first():
    list(parse_fasta([">fasta1\n", ">fasta2\n", "AGTC\n"]))

@nose.tools.raises(FASTAError)
def test_parse_fasta__empty_record__middle():
    lines = [">fasta0\n", "ACGT\n", ">fasta1\n", ">fasta2\n", "AGTC\n"]
    list(parse_fasta(lines))

@nose.tools.raises(FASTAError)
def test_parse_empty_record_last():
    lines = [">fasta1\n", "ACGT\n", ">fasta2\n"]
    list(parse_fasta(lines))

@nose.tools.raises(FASTAError)
def test_parse_fasta__missing_name__alone():
    lines = ["ACGT\n"]
    list(parse_fasta(lines))

@nose.tools.raises(FASTAError)
def test_parse_fasta__missing_name__with_others():
    lines = ["ACGT\n", ">Foo\n", "ACGGTA\n"]
    list(parse_fasta(lines))

@nose.tools.raises(FASTAError)
def test_parse_fasta__empty_name__alone():
    lines = [">\n", "ACGT\n"]
    list(parse_fasta(lines))

@nose.tools.raises(FASTAError)
def test_parse_fasta__empty_name__with_others():
    lines = [">\n", "ACGT\n", ">Foo\n", "ACGGTA\n"]
    list(parse_fasta(lines))




################################################################################
################################################################################
## Tests for 'read_fasta'

def test_read_fasta__uncompressed():
    expected = [(("This_is_FASTA!", None), "ACGTN"),
                (("This_is_ALSO_FASTA!", None), "CGTNA")]
    results  = list(read_fasta("tests/data/fasta_file.fasta"))
    assert_equals(results, expected)

def test_read_fasta__compressed_gz():
    expected = [(("This_is_GZipped_FASTA!", None), "ACGTN"),
                (("This_is_ALSO_GZipped_FASTA!", None), "CGTNA")]
    results  = list(read_fasta("tests/data/fasta_file.fasta.gz"))
    assert_equals(results, expected)

def test_read_fasta__compressed_bz2():
    expected = [(("This_is_BZ_FASTA!", None), "CGTNA"),
                (("This_is_ALSO_BZ_FASTA!", None), "ACGTN")]
    results  = list(read_fasta("tests/data/fasta_file.fasta.bz2"))
    assert_equals(results, expected)
