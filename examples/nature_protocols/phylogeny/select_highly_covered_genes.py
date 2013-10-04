#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import sys

from pypeline.common.formats.msa import \
     MSA

# Minimum fraction of bases that must be covered across all sequences
_FRAC_BASES_COVERED = 0.8
# Exclude the following sequence(s) when calculating covered bases
_EXCLUDE_SAMPLES = frozenset(("Pi_nucl",))


def _is_sufficently_covered(filepath):
    msa = MSA.from_file(filepath)
    if msa.seqlen() % 3:
        return False

    filtered_msa = msa.exclude(_EXCLUDE_SAMPLES)

    total_bases_not_covered = 0
    for fasta_record in filtered_msa:
        total_bases_not_covered += fasta_record.sequence.upper().count("N")
        total_bases_not_covered += fasta_record.sequence.count("-")

    frac_covered = (1.0 - total_bases_not_covered / float(len(filtered_msa) * filtered_msa.seqlen()))
    return frac_covered >= _FRAC_BASES_COVERED


def main(argv):
    for root_dir in argv:
        for filename in os.listdir(root_dir):
            if filename.endswith(".afa"):
                if _is_sufficently_covered(os.path.join(root_dir, filename)):
                    sequence_name, _ = filename.rsplit(".", 1)
                    print(sequence_name)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
