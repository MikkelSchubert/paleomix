import gzip
from collections import defaultdict

from pypeline.common.sequences import split
from pypeline.common.formats.fasta import read_fasta, parse_fasta


class MSAError(RuntimeError):
    pass


def split_msa(msa, split_by = "123"):
    _validate_msa_lengths(msa)
    if not split_by:
        raise TypeError("No partitions to split by specified")


    results = {}
    for key in split_by:
        results[key] = dict((name, {}) for name in msa)

    for (name, sequence) in msa.iteritems():
        for (key, partition) in split(sequence, split_by).iteritems():
            results[key][name] = partition

    return results
    
    

def join_msa(*msas):
    if not msas:
        raise TypeError("join_msa takes at least 1 argument (0 given)")

    results = defaultdict(list)
    keywords = set(msas[0])
    for msa in msas:
        _validate_msa_lengths(msa)

        if set(msa) != keywords:
            raise MSAError("Cannot join MSAs with mismatching sequences")
        
        for (name, sequence) in msa.iteritems():
            results[name].append(sequence)

    return dict((key, "".join(value)) for (key, value) in results.iteritems())
    

def parse_msa(lines):
    msa = {}
    for (name, sequence) in parse_fasta(lines):
        if name in msa:
            raise MSAError("Duplicate names found, cannot be represented as MSA")
        msa[name] = sequence


    return _validate_msa_lengths(msa)


def read_msa(filename):
    func = gzip.open if filename.endswith(".gz") else open
    fasta_file = func(filename, "r")
    try:
        return parse_msa(iter(fasta_file))
    finally:
        fasta_file.close()


def _validate_msa_lengths(msa):
    if len(set(len(seq) for seq in msa.itervalues())) > 1:
        raise MSAError("Sequence lengths differ in MSA")

    return msa
