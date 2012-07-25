import gzip
import types
from collections import defaultdict

from pypeline.common.sequences import split
from pypeline.common.formats.fasta import parse_fasta


class MSAError(RuntimeError):
    pass


def split_msa(msa, split_by = "123"):
    validate_msa(msa)
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
    validate_msa(*msas)
    results = defaultdict(list)
    for msa in msas:
        for (name, sequence) in msa.iteritems():
            results[name].append(sequence)

    return dict((key, "".join(value)) for (key, value) in results.iteritems())
    

def parse_msa(lines):
    msa = {}
    for (name, sequence) in parse_fasta(lines):
        if name in msa:
            raise MSAError("Duplicate names found, cannot be represented as MSA")
        msa[name] = sequence

    validate_msa(msa)
    return msa


def read_msa(filename):
    func = gzip.open if filename.endswith(".gz") else open
    fasta_file = func(filename, "r")
    try:
        return parse_msa(iter(fasta_file))
    finally:
        fasta_file.close()


def validate_msa(*msas):
    if not msas:
        raise TypeError("No MSAs given as arguments")

    keywords = set(msas[0])
    for msa in msas:
        if not msa:
            raise MSAError("MSA with no sequences found")
        elif not all((name and isinstance(name, types.StringTypes)) for name in msa):
            raise MSAError("Names in MSA must be non-empty strings")
        elif len(set(len(seq) for seq in msa.itervalues())) != 1:
            raise MSAError("MSA contains sequences of differing lengths")
        elif set(msa) != keywords:
            raise MSAError("MSAs contain mismatching sequences")
