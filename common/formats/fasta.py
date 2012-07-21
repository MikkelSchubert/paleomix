from __future__ import print_function

import sys
import gzip

from common.utilities import grouper, split_before




def wrap_fasta(name, sequence):
    """Process a printable FASTA sequence, wrapping long sequences at 60 chars."""
    lines = [">%s" % (name,)]
    for line in grouper(60, sequence, fillvalue = ""):
        lines.append("".join(line))
    lines.append("")
    return "\n".join(lines)


def print_fasta(name, sequence, file = sys.stdout):
    """Prints a FASTA sequence (iterable), wrapping long sequences at 60 chars."""
    print(wrap_fasta(name, sequence), file = file, end = "")


def parse_fasta(lines):
    """Parses FASTA sequences found in a sequence of lines, and returns
    a tuple containing the name and full sequence of each FASTA record.
    No assumptions are made about the line-lengths."""
    lines = (line.rstrip() for line in lines)
    for fragment in split_before(lines, lambda v: v.startswith(">")):
        name = fragment[0]
        if not name.startswith(">"):
            raise ValueError("Unnamed FASTA record")
        elif len(fragment) == 1:
            raise ValueError("FASTA record does not contain sequence.")

        yield (name[1:], "".join(fragment[1:]))


def read_fasta(filename):
    """Reads an unindexed FASTA file, returning a sequence of 
    tuples containing the name and sequence of each entry in
    the file. The FASTA file may be GZIP compressed."""
    
    func = gzip.open if filename.endswith(".gz") else open
    fasta_file = func(filename, "r")
    try:
        for record in parse_fasta(iter(fasta_file)):
            yield record
    finally:
        fasta_file.close()
