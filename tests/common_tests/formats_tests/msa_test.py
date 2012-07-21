import nose.tools

from pypeline.common.formats.msa import *




################################################################################
################################################################################
## Tests for 'parse_msa'

def test_parse_msa__unique_names():
    assert False

@nose.tools.raises(MSAError)
def test_parse_msa__duplicate_names():
    pass

@nose.tools.raises(MSAError)
def test_parse_msa__mismatched_lengths():
    pass




################################################################################
################################################################################
## Tests for 'read_msa'

def test_read_msa__uncompressed():
    expected = {"This is FASTA!" : "ACGTN",
                "This is ALSO FASTA!" : "CGTNA"}
    results  = read_msa("tests/data/fasta_file.fasta")

    assert results == expected

def test_read_msa__compressed():
    expected = {"This is GZipped FASTA!" : "ACGTN",
                "This is ALSO GZipped FASTA!" : "CGTNA"}
    results  = read_msa("tests/data/fasta_file.fasta.gz")

    assert results == expected
