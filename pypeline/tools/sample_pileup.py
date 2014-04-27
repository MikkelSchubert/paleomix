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
from __future__ import print_function
from __future__ import with_statement

import sys
import random
import itertools
import argparse

import pysam

import pypeline.common.text as text
import pypeline.common.sequences as sequences
from pypeline.common.formats.fasta import FASTA


class Pileup:
    def __init__(self, line):
        fields = line.split("\t")

        self.position = int(fields[1]) - 1
        self.reference = fields[2]
        self.depth = int(fields[3])
        self.observed = fields[4]
        self.qualities = fields[5]


def _untrusted_positions(pileups, distance):
    """Returns a set of positions that are either directly covered by, or
    adjacent to indels, given some arbitrary distance. These positions are
    considered 'untrusted', and bases that fall within those positions are
    ignored.
    """
    positions = set()
    for pileup in pileups:
        for (index, current) in enumerate(pileup.observed):
            if current == "+":
                length = 0
            elif current == "-":
                length = int(pileup.observed[index + 1])
            else:
                continue

            # Inclusive start/end positions for blacklisted bases
            start = pileup.position - distance
            end = pileup.position + distance + length

            positions.update(xrange(start, end + 1))

    return positions


def sample_position(pileup):
    skip = 0
    bases = []
    observed = pileup.observed.upper()
    for (current, lookahead) in zip(observed, observed[1:] + "N"):
        if skip:
            skip -= 1
        elif current in "+-":
            skip = int(lookahead) + 1
        elif current in ".,":
            bases.append(pileup.reference)
        elif current in "ACGT":
            bases.append(current)
        elif current == "^":
            skip = 1
        elif current not in "$*N":
            assert False, current

    if not bases:
        return "N"

    return random.choice(bases)


def build_region(options, genotype, bed):
    # Note that bed.end is a past-the-end coordinate
    start = max(0, bed.start - options.padding)
    end = bed.end + options.padding

    pileups = []
    if bed.contig in genotype.contigs:
        # fetch raises a ValueError if the VCF does not contain any entries for
        # the specified contig, which occur due to low coverage.
        pileups = [Pileup(pileup)
                   for pileup in genotype.fetch(bed.contig, start, end)]
    untrusted = _untrusted_positions(pileups, options.min_distance_to_indels)

    sequence = ["N"] * (end - start)
    for pileup in pileups:
        if pileup.position not in untrusted:
            sequence[pileup.position - start] = sample_position(pileup)

    offset = bed.start - start
    length = bed.end - bed.start
    return sequence[offset: offset + length]


def build_genes(options, genotype, regions):
    def keyfunc(bed):
        return (bed.contig, bed.name, bed.start)
    regions.sort(key=keyfunc)

    for (gene, beds) in itertools.groupby(regions, lambda x: x.name):
        sequence, beds = [], tuple(beds)
        for bed in beds:
            sequence.extend(build_region(options, genotype, bed))
        sequence = "".join(sequence)

        if any((bed.strand == "-") for bed in beds):
            assert all((bed.strand == "-") for bed in beds)

            sequence = sequences.reverse_complement(sequence)

        yield (gene, sequence)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--genotype", help="Tabix indexed pileup file.",
                        required=True)
    parser.add_argument("--intervals", help="BED file.", required=True)
    parser.add_argument("--padding", type=int, default=10,
                        help="Number of bases to expand intervals, when "
                             "filtering based on adjacent indels [%default]")
    parser.add_argument("--min-distance-to-indels", type=int, default=5,
                        help="Variants closer than this distance from indels "
                             "are filtered [%default].")
    args = parser.parse_args(argv)

    genotype = pysam.Tabixfile(args.genotype)
    with open(args.intervals) as bed_file:
        intervals = text.parse_lines_by_contig(bed_file, pysam.asBed())

    for (_, beds) in sorted(intervals.items()):
        for (name, sequence) in build_genes(args, genotype, beds):
            FASTA(name, None, sequence).write(sys.stdout)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
