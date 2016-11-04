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
from __future__ import print_function
from __future__ import with_statement

import argparse
import collections
import heapq
import itertools
import random
import sys

import pysam

import paleomix.common.sequences as sequences
import paleomix.common.text as text

from paleomix.common.bedtools import BEDRecord
from paleomix.common.formats.fasta import FASTA


# Regions are genotyped in chunks of 1 MB strings; this is done to reduce
# overhead, as storing 1 MB of chars individually in a list adds about 45 MB of
# overhead. A 100MB chromosome would therefore require ~4.5GB.
_CHUNK_SIZE = 2 ** 20


class Pileup(object):
    def __init__(self, line):
        fields = line.split("\t")

        self.position = int(fields[1]) - 1
        self.reference = fields[2]
        self.depth = int(fields[3])
        self.observed = fields[4]


class PileupRegion(object):
    def __init__(self, tabix, chrom, start, end, padding=0,
                 dist_to_indels=0):
        assert padding >= 0, "Padding must be >= 0, not %i" % (padding,)
        self._tabix = tabix
        self._chrom = chrom
        self._start = start
        self._end = end
        self._padding = padding
        self._distance_to_indels = dist_to_indels

    def __iter__(self):
        if self._distance_to_indels <= 0:
            for pileup in self._tabix.fetch(reference=self._chrom,
                                            start=self._start,
                                            end=self._end):
                yield (pileup.position, self._sample_position(pileup))

        # Note that bed.end is a past-the-end coordinate
        start = max(0, self._start - self._padding)
        end = self._end + self._padding

        region = self._tabix.fetch(reference=self._chrom,
                                   start=start,
                                   end=end)
        pileup_buffer = collections.deque()
        blacklist = []

        # Fill buffer, and detect blacklisted sites due to indels
        for _ in xrange(max(self._padding, self._distance_to_indels) * 2):
            self._add_to_buffer(region, pileup_buffer, blacklist)

        while pileup_buffer:
            position, nucleotide = pileup_buffer.popleft()
            while blacklist and blacklist[0] < position:
                heapq.heappop(blacklist)

            if not blacklist or blacklist[0] != position:
                if self._start <= position < self._end:
                    yield (position, nucleotide)

            self._add_to_buffer(region, pileup_buffer, blacklist)

    def _add_to_buffer(self, region, pileup_buffer, blacklist):
        try:
            pileup = Pileup(region.next())
            self._collect_indels(pileup, blacklist)
            pileup_buffer.append((pileup.position,
                                  self._sample_position(pileup)))
        except StopIteration:
            pass

    def _collect_indels(self, pileup, blacklist):
        previous = None
        for (index, current) in enumerate(pileup.observed):
            if previous == '^':
                previous = current
                continue

            previous = current
            if current == "+":
                # Insertions do not themselves cover any bases
                length = 0
            elif current == "-":
                len_slice = itertools.islice(pileup.observed, index + 1, None)
                digits = "".join(itertools.takewhile(str.isdigit, len_slice))
                length = int(digits)
            else:
                continue

            # Distance is defined as sites overlapping INDELs having distance
            # 0, sites adjacent to INDELS have distance 1, etc. Note that the
            # INDEL starts on the next position of the current row.
            start = pileup.position - self._distance_to_indels + 2
            end = pileup.position + self._distance_to_indels + length
            for position in xrange(start, end):
                heapq.heappush(blacklist, position)

    @classmethod
    def _sample_position(cls, pileup):
        skip = 0
        bases = []
        observed = pileup.observed.upper()
        for current in observed:
            if skip:
                skip -= 1
            elif current in ".,":
                bases.append(pileup.reference)
            elif current in "ACGT":
                bases.append(current)
            elif current == "^":
                skip = 1
            elif current not in "$*N+-1234567890":
                assert False, current

        if not bases:
            return "N"

        return random.choice(bases)


def build_region(options, genotype, bed):
    # 'fetch' raises a ValueError if the VCF does not contain any entries for
    # the specified contig, which can occur due to low coverage.
    if bed.contig in genotype.contigs:
        region = PileupRegion(genotype,
                              chrom=bed.contig,
                              start=bed.start,
                              end=bed.end,
                              padding=options.padding,
                              dist_to_indels=options.min_distance_to_indels)

        remaining_length = (bed.end - bed.start)
        offset = bed.start

        # Genotyping is done in chunks, so that these can be reduced to strings
        # and thereby reduce the memory requirements for larger regions.
        chunk = ["N"] * min(_CHUNK_SIZE, remaining_length)

        for position, nucleotide in region:
            while position >= offset + len(chunk):
                yield "".join(chunk)
                remaining_length -= len(chunk)
                offset += len(chunk)
                chunk = ["N"] * min(_CHUNK_SIZE, remaining_length)

            chunk[position - offset] = nucleotide

        while offset < bed.end:
            yield "".join(chunk)
            remaining_length -= len(chunk)
            offset += len(chunk)
            chunk = ["N"] * min(_CHUNK_SIZE, remaining_length)

        yield "".join(chunk)
    else:
        yield "N" * (bed.end - bed.start)


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
    prog = "paleomix sample_pileup"
    usage = "%s [options] --genotype in.vcf --intervals in.bed > out.fasta" \
        % (prog,)

    parser = argparse.ArgumentParser(prog=prog, usage=usage)
    parser.add_argument("--genotype", help="Tabix indexed pileup file.",
                        required=True, metavar="PILEUP")
    parser.add_argument("--intervals", help="BED file.", required=True,
                        metavar="BED")
    parser.add_argument("--padding", type=int, default=10, metavar="BED",
                        help="Number of bases to expand intervals, when "
                             "filtering based on adjacent indels "
                             "[%(default)s]")
    parser.add_argument("--min-distance-to-indels", type=int, default=5,
                        help="Variants closer than this distance from indels "
                             "are filtered; set to a negative value to "
                             "disable [%(default)s].")
    args = parser.parse_args(argv)

    genotype = pysam.Tabixfile(args.genotype)
    with open(args.intervals) as bed_file:
        intervals = text.parse_lines_by_contig(bed_file, BEDRecord)

    for (_, beds) in sorted(intervals.items()):
        for (name, sequence) in build_genes(args, genotype, beds):
            FASTA(name, None, sequence).write(sys.stdout)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
