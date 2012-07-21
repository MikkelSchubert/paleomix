import gzip
from collections import defaultdict

from pypeline.common.formats.fasta import read_fasta, parse_fasta


class MSAError(RuntimeError):
    pass



def join_msa(*msas):
    if not msas:
        raise TypeError("join_msa takes at least 1 argument (0 given)")

    results = defaultdict(list)
    keywords = set(msas[0])
    for msa in msas:
        if set(msa) != keywords:
            raise ValueError("Cannot join MSAs with mismatching sequences")
        
        for (name, sequence) in msa.iteritems():
            results[name].append(sequence)

    return dict((key, "".join(value)) for (key, value) in results.iteritems())
    

def parse_msa(lines):
    msa = dict(parse_fasta(lines))
    if len(set(len(seq) for seq in msa.itervalues())) != 1:
        raise MSAError("Sequence lengths differ in MSA")

    return msa


def read_msa(filename):
    func = gzip.open if filename.endswith(".gz") else open
    fasta_file = func(filename, "r")
    try:
        return parse_msa(iter(fasta_file))
    finally:
        fasta_file.close()
