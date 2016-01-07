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
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
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
import optparse

from paleomix.common.formats.msa import \
     MSA


def _is_sufficently_covered(filepath, min_coverage):
    msa = MSA.from_file(filepath)
    if msa.seqlen() % 3:
        return False

    total_bases_not_covered = 0
    for fasta_record in msa:
        total_bases_not_covered += fasta_record.sequence.upper().count("N")
        total_bases_not_covered += fasta_record.sequence.count("-")

    total_bases = float(len(msa) * msa.seqlen())
    frac_covered = 1.0 - total_bases_not_covered / total_bases
    return frac_covered >= min_coverage


def main(argv):
    usage = "%prog [options] <PATH>"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--min-coverage", default=0.8, type=float,
                      help="Minimum fraction of called bases in a MSA "
                           "[%default]")

    config, args = parser.parse_args(argv)
    if not args:
        parser.print_usage()
        return 1

    for root_dir in args:
        for filename in os.listdir(root_dir):
            if filename.endswith(".afa"):
                fpath = os.path.join(root_dir, filename)
                if _is_sufficently_covered(fpath, config.min_coverage):
                    sequence_name, _ = filename.rsplit(".", 1)
                    print(sequence_name)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
