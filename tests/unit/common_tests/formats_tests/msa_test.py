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
import nose.tools
from nose.tools import assert_equal

from pypeline.common.formats.msa import *



################################################################################
################################################################################
## Tests for 'join_msa'

_JOIN_MSA_1 = { "nc" : "ACG", "nm" : "TGA", "miRNA" : "UCA" }
_JOIN_MSA_2 = { "nc" : "TGA", "nm" : "CTT", "miRNA" : "GAC" }
_JOIN_MSA_3 = { "nc" : "AAG", "nm" : "GAG", "miRNA" : "CAU" }

def test_join_msa__single_msa():
    result = join_msa(_JOIN_MSA_1)
    assert_equal(result, _JOIN_MSA_1)

def test_join_msa__two_msa():
    expected = { "nc" : "ACGTGA", "nm" : "TGACTT", "miRNA" : "UCAGAC" }
    result = join_msa(_JOIN_MSA_1, _JOIN_MSA_2)
    assert_equal(result, expected)

def test_join_msa__three_msa():
    expected = { "nc" : "ACGTGAAAG", 
                 "nm" : "TGACTTGAG", 
                 "miRNA" : "UCAGACCAU" }
    result = join_msa(_JOIN_MSA_1, _JOIN_MSA_2, _JOIN_MSA_3)
    assert_equal(result, expected)

@nose.tools.raises(TypeError)
def test_join_msa__missing_arguments():
    join_msa()

@nose.tools.raises(MSAError)
def test_join_msa__missing_names_first():
    msa_1 = dict(_JOIN_MSA_1)
    msa_2 = dict(_JOIN_MSA_2)
    msa_1.pop(msa_1.keys()[0])
    join_msa(msa_1, msa_2)

@nose.tools.raises(MSAError)
def test_join_msa__missing_names_second():
    msa_1 = dict(_JOIN_MSA_1)
    msa_2 = dict(_JOIN_MSA_2)
    msa_2.pop(msa_2.keys()[0])
    join_msa(msa_1, msa_2)

@nose.tools.raises(MSAError)
def test_join_msa__differing_lengths():
    msa_1 = dict(_JOIN_MSA_1)
    msa_2 = dict(_JOIN_MSA_2)
    msa_1["nc"] = "AC"
    join_msa(msa_1, msa_2)

@nose.tools.raises(MSAError)
def test_join_msa__empty_name():
    join_msa({"" : "A"}, {"" : "T"})

@nose.tools.raises(MSAError)
def test_join_msa__non_string_name():
    join_msa({1 : "A"}, {1 : "T"})



################################################################################
################################################################################
## Tests for 'parse_msa'

def test_parse_msa__single_entry():
    lines = [">seq1", "ACG"]
    assert_equal(parse_msa(lines), {"seq1" : "ACG"})

def test_parse_msa__two_entries():
    lines = [">seq1", "ACG", ">seq2", "TGA"]
    assert_equal(parse_msa(lines), {"seq1" : "ACG", "seq2" : "TGA"})

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
    expected = {"This is FASTA!" : "ACGTN",
                "This is ALSO FASTA!" : "CGTNA"}
    results  = read_msa("tests/data/fasta_file.fasta")
    assert_equal(results, expected)

def test_read_msa__compressed():
    expected = {"This is GZipped FASTA!" : "ACGTN",
                "This is ALSO GZipped FASTA!" : "CGTNA"}
    results  = read_msa("tests/data/fasta_file.fasta.gz")
    assert_equal(results, expected)


################################################################################
################################################################################
## Tests for 'split_msa'

def test_split_msa__single_group():
    msa = {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}
    expected = {'1' : {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}}
    assert_equal(split_msa(msa, "111"), expected)

def test_split_msa__two_groups(): 
    msa = {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}
    expected = {"1" : {"seq1" : "ACCA", "seq2" : "GATG"},
                "2" : {"seq1" : "GT",   "seq2" : "GA"}}
    assert_equal(split_msa(msa, "112"), expected)

def test_split__three_groups():
    msa = {"seq1" : "ACGCAT", "seq2" : "GAGTGA"}
    expected = {"1" : {"seq1" : "AC", "seq2" : "GT"},
                "2" : {"seq1" : "CA", "seq2" : "AG"},
                "3" : {"seq1" : "GT", "seq2" : "GA"}}
    assert_equal(split_msa(msa, "123"), expected)
    assert_equal(split_msa(msa), expected)

def test_split__empty_group():
    msa = {"seq1" : "AC", "seq2" : "GA"}
    expected = {"1" : {"seq1" : "A", "seq2" : "G"},
                "2" : {"seq1" : "C", "seq2" : "A"},
                "3" : {"seq1" : "",  "seq2" : ""}}
    assert_equal(split_msa(msa), expected)

def test_split__partial_group():
    msa = {"seq1" : "ACGCA", "seq2" : "GAGTG"}
    expected = {"1" : {"seq1" : "AC", "seq2" : "GT"},
                "2" : {"seq1" : "CA", "seq2" : "AG"},
                "3" : {"seq1" : "G", "seq2" : "G"}}
    assert_equal(split_msa(msa), expected)


@nose.tools.raises(MSAError)
def test_split_msa__empty_msa():
    split_msa({})

@nose.tools.raises(TypeError)
def test_split_msa__no_split_by():
    split_msa({"seq1" : "ACG", "seq2" : "GAT"}, split_by = "")

@nose.tools.raises(MSAError)
def test_split_msa__mismatching_lengths():
    split_msa({"seq1" : "ACG", "seq2" : "GA"})

@nose.tools.raises(MSAError)
def test_split_msa__empty_name():
    split_msa({"" : "ACT"}, {"" : "TGA"})

@nose.tools.raises(MSAError)
def test_split_msa__non_string_name():
    split_msa({1 : "ACT"}, {1 : "TGA"})
