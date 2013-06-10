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
import os
import StringIO

import nose.tools
from nose.tools import assert_equal # pylint: disable=E0611

from tests.common.utils import require_call, with_temp_folder
from pypeline.common.formats.msa import \
     split_msa, \
     join_msa, \
     parse_msa, \
     read_msa, \
     print_msa, \
     write_msa, \
     validate_msa, \
     MSAError

_VALIDATION_PATH = "pypeline.common.formats.msa.validate_msa"


################################################################################
################################################################################
## Tests for 'join_msa'

_JOIN_MSA_1 = { "nc" : "ACG", "nm" : "TGA", "miRNA" : "UCA" }
_JOIN_MSA_2 = { "nc" : "TGA", "nm" : "CTT", "miRNA" : "GAC" }
_JOIN_MSA_3 = { "nc" : "AAG", "nm" : "GAG", "miRNA" : "CAU" }

def test_join_msa__single_msa():
    with require_call(_VALIDATION_PATH, args = [_JOIN_MSA_1]):
        result = join_msa(_JOIN_MSA_1)
        assert_equal(result, _JOIN_MSA_1)

def test_join_msa__two_msa():
    expected = { "nc" : "ACGTGA", "nm" : "TGACTT", "miRNA" : "UCAGAC" }
    with require_call(_VALIDATION_PATH, args = [_JOIN_MSA_1, _JOIN_MSA_2]):
        result = join_msa(_JOIN_MSA_1, _JOIN_MSA_2)
        assert_equal(result, expected)

def test_join_msa__three_msa():
    expected = { "nc" : "ACGTGAAAG",
                 "nm" : "TGACTTGAG",
                 "miRNA" : "UCAGACCAU" }
    with require_call(_VALIDATION_PATH, args = [_JOIN_MSA_1, _JOIN_MSA_2, _JOIN_MSA_3]):
        result = join_msa(_JOIN_MSA_1, _JOIN_MSA_2, _JOIN_MSA_3)
        assert_equal(result, expected)

@nose.tools.raises(TypeError)
def test_join_msa__missing_arguments():
    join_msa()




################################################################################
################################################################################
## Tests for 'parse_msa'

def test_parse_msa__single_entry():
    lines  = [">seq1", "ACG"]
    result = {"seq1" : "ACG"}
    with require_call(_VALIDATION_PATH, args = [result]):
        assert_equal(parse_msa(lines), result)

def test_parse_msa__single_entry_with_meta():
    lines  = [">seq1 Meta info", "ACG"]
    result_msa  = {"seq1" : "ACG"}
    result_meta = {"seq1" : "Meta info"}
    with require_call(_VALIDATION_PATH, args = [result_msa]):
        msa, meta = parse_msa(lines, read_meta = True)
        assert_equal(msa, result_msa)
        assert_equal(meta, result_meta)

def test_parse_msa__two_entries():
    lines  = [">seq1", "ACG", ">seq2", "TGA"]
    result = {"seq1" : "ACG", "seq2" : "TGA"}
    with require_call(_VALIDATION_PATH, args = [result]):
        assert_equal(parse_msa(lines), result)

def test_parse_msa__two_entries_with_meta():
    lines  = [">seq1", "ACG", ">seq2 Second meta", "TGA"]
    result_msa = {"seq1" : "ACG", "seq2" : "TGA"}
    result_meta = {"seq1" : None, "seq2" : "Second meta"}
    with require_call(_VALIDATION_PATH, args = [result_msa]):
        msa, meta = parse_msa(lines, read_meta = True)
        assert_equal(msa, result_msa)
        assert_equal(meta, result_meta)

@nose.tools.raises(MSAError)
def test_parse_msa__duplicate_names():
    parse_msa([">seq1", "ACG", ">seq1", "TGA"])

@nose.tools.raises(MSAError)
def test_parse_msa__mismatched_lengths():
    parse_msa([">seq1", "ACG", ">seq1", "TGAN"])

@nose.tools.raises(MSAError)
def test_parse_msa__empty_name():
    parse_msa([">", "ACG", ">seq1", "TGAN"])




################################################################################
################################################################################
## Tests for 'read_msa'

def test_read_msa__uncompressed():
    expected = {"This_is_FASTA!" : "ACGTN",
                "This_is_ALSO_FASTA!" : "CGTNA"}
    with require_call(_VALIDATION_PATH, args = [expected]):
        results  = read_msa("tests/data/fasta_file.fasta")
        assert_equal(results, expected)

def test_read_msa__compressed_gz():
    expected = {"This_is_GZipped_FASTA!" : "ACGTN",
                "This_is_ALSO_GZipped_FASTA!" : "CGTNA"}
    with require_call(_VALIDATION_PATH, args = [expected]):
        results  = read_msa("tests/data/fasta_file.fasta.gz")
        assert_equal(results, expected)

def test_read_msa__compressed_bz2():
    expected = {"This_is_BZ_FASTA!" : "CGTNA",
                "This_is_ALSO_BZ_FASTA!" : "ACGTN"}
    with require_call(_VALIDATION_PATH, args = [expected]):
        results  = read_msa("tests/data/fasta_file.fasta.bz2")
        assert_equal(results, expected)




################################################################################
################################################################################
## Tests for 'split_msa'

def test_split_msa__single_group():
    msa = {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}
    expected = {'1' : {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}}
    with require_call(_VALIDATION_PATH, args = [msa]):
        assert_equal(split_msa(msa, "111"), expected)

def test_split_msa__two_groups():
    msa = {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}
    expected = {"1" : {"seq1" : "ACCA", "seq2" : "GATG"},
                "2" : {"seq1" : "GT",   "seq2" : "GA"}}
    with require_call(_VALIDATION_PATH, args = [msa]):
        assert_equal(split_msa(msa, "112"), expected)

def test_split__three_groups():
    msa = {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}
    expected = {"1" : {"seq1" : "AC", "seq2" : "GT"},
                "2" : {"seq1" : "CA", "seq2" : "AG"},
                "3" : {"seq1" : "GT", "seq2" : "GA"}}
    with require_call(_VALIDATION_PATH, args = [msa]):
        assert_equal(split_msa(msa, "123"), expected)
    with require_call(_VALIDATION_PATH, args = [msa]):
        assert_equal(split_msa(msa), expected)

def test_split__empty_group():
    msa = {"seq1" : "AC", "seq2" : "GA"}
    expected = {"1" : {"seq1" : "A", "seq2" : "G"},
                "2" : {"seq1" : "C", "seq2" : "A"},
                "3" : {"seq1" : "",  "seq2" : ""}}
    with require_call(_VALIDATION_PATH, args = [msa]):
        assert_equal(split_msa(msa), expected)

def test_split__partial_group():
    msa = {"seq1" : "ACGCA", "seq2" : "GAGTG"}
    expected = {"1" : {"seq1" : "AC", "seq2" : "GT"},
                "2" : {"seq1" : "CA", "seq2" : "AG"},
                "3" : {"seq1" : "G", "seq2" : "G"}}
    with require_call(_VALIDATION_PATH, args = [msa]):
        assert_equal(split_msa(msa), expected)


@nose.tools.raises(MSAError)
def test_split_msa__empty_msa():
    split_msa({})

@nose.tools.raises(TypeError)
def test_split_msa__no_split_by():
    split_msa({"seq1" : "ACG", "seq2" : "GAT"}, split_by = "")




################################################################################
################################################################################
## Tests for 'write_msa'

@with_temp_folder
def test_write_msa(temp_folder):
    msa = {"seq1" : "ACGTA", "seq2" : "CGTAC"}
    fname = os.path.join(temp_folder, "out.afa")
    with require_call(_VALIDATION_PATH, args = [msa]):
        write_msa(msa, fname)
    with require_call(_VALIDATION_PATH, args = [msa]):
        assert_equal(read_msa(fname), msa)




################################################################################
################################################################################
## Tests for 'print_msa'

def test_print_fasta__complete_line_test():
    msa       = {"barfoo" : "ACGATA" * 10 + "CGATAG" * 5,
                 "foobar" : "CGAATG" * 10 + "TGTCAT" * 5}
    expected  = ">barfoo\n%s\n%s\n" % ("ACGATA" * 10, "CGATAG" * 5)
    expected += ">foobar\n%s\n%s\n" % ("CGAATG" * 10, "TGTCAT" * 5)
    stringf = StringIO.StringIO()
    with require_call(_VALIDATION_PATH, args = [msa]):
        print_msa(msa, stringf)
    assert_equal(stringf.getvalue(), expected)



################################################################################
################################################################################
## Tests for 'validate_msa'


@nose.tools.raises(MSAError)
def test_validate_msa__missing_names_first():
    msa_1 = dict(_JOIN_MSA_1)
    msa_2 = dict(_JOIN_MSA_2)
    msa_1.pop(msa_1.keys()[0])
    validate_msa(msa_1, msa_2)

@nose.tools.raises(MSAError)
def test_validate_msa__missing_names_second():
    msa_1 = dict(_JOIN_MSA_1)
    msa_2 = dict(_JOIN_MSA_2)
    msa_2.pop(msa_2.keys()[0])
    validate_msa(msa_1, msa_2)

@nose.tools.raises(MSAError)
def test_validate_msa__differing_lengths():
    msa_1 = dict(_JOIN_MSA_1)
    msa_2 = dict(_JOIN_MSA_2)
    msa_1["nc"] = "AC"
    validate_msa(msa_1, msa_2)

@nose.tools.raises(MSAError)
def test_validate_msa__empty_name():
    validate_msa({"" : "A"}, {"" : "T"})

@nose.tools.raises(MSAError)
def test_validate_msa__non_string_name():
    validate_msa({1 : "A"}, {1 : "T"})
