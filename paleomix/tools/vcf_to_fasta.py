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
"""Creats consensus sequence from a VCF and BED file:
vcf_to_fasta --intervals PATH_TO.bed --genotype PATH_TO.vcf

The VCF file is expected to have been created using SAMTools, and makes use
of SAMTools specific fields. For each region in the BED file, the script will
create a consensus sequence containing all basepairs which are filtered as
'PASS' or '.' in the VCF file; at each site, the most likely phenotype is
selected using the PL fields, and heterozygous sites encoded using UIPAC codes.

BED features with the same name are merged into single sequences. To ensure
that indels are called near sequence termini, the script expects the VCF file
to contain a certain amount of padding around the regions of interest.

"""
from __future__ import print_function

import argparse
import copy
import itertools
import os
import sys
import re

import pysam

import paleomix.common.vcfwrap as vcfwrap
import paleomix.common.text as text
import paleomix.common.sequences as sequences
import paleomix.common.utilities as utilities

from paleomix.common.bedtools import BEDRecord


# Max number of positions to keep in memory / genotype at once
_SEQUENCE_CHUNK = 1024 * 1024  # 1kbp
# Number of columns per line in FASTA sequences
_FASTA_COLUMNS = 60

_VCF_DICT = re.compile("##(.*)=<(.*)>")


###############################################################################
###############################################################################
# Utility functions

def flush_fasta(sequence):
    """Takes a FASTA sequence as a string, fragments it into lines of exactly
    _FASTA_COLUMNS chars (e.g. 60), and prints all complete lines. The final
    incomplete line (if any) is returned.

    """
    for seq_frag in utilities.fragment(_FASTA_COLUMNS, sequence):
        if len(seq_frag) < _FASTA_COLUMNS:
            return seq_frag
        print(seq_frag)
    return ""


def split_beds(beds, size=_SEQUENCE_CHUNK):
    """Takes a list of beds, and splits each bed into chunks that are at most
    'size' bp long. The resulting (smaller) beds are returned as a new list.

    """
    results = []
    for bed in beds:
        for start in xrange(bed.start, bed.end, size):
            end = min(start + size, bed.end)
            split_bed = copy.copy(bed)
            split_bed.start = start
            split_bed.end = end
            results.append(split_bed)
    return results


###############################################################################
###############################################################################
# Genotyping functions

def add_snp(options, snp, position, sequence):
    if snp.alt != ".":
        genotype = "".join(vcfwrap.get_ml_genotype(snp, options.nth_sample))
        encoded = sequences.encode_genotype(genotype)
    else:
        encoded = snp.ref
    sequence[position] = encoded


def add_indel(options, bed, indel, sequence):
    if indel.alt == ".":
        return

    genotype = vcfwrap.get_ml_genotype(indel, options.nth_sample)
    if genotype[0] != genotype[1]:
        # No way to represent heterozygous indels
        return
    elif genotype[0] == "N":
        # No most likely genotype
        return

    # Note that bed.end is a past-the-end coordinate
    start = max(0, bed.start - options.padding)

    # FIXME: parse_indel only supports a single 'alt' values
    indel.alt = genotype[0]
    indel = vcfwrap.parse_indel(indel)
    if indel.in_reference:
        del_start = max(indel.pos + 1, bed.start)
        del_end = min(indel.pos + 1 + len(indel.what), bed.end)

        if del_start >= del_end:
            # Deletion does not cover any bases of interest
            return
        elif options.whole_codon_indels_only:
            if (del_end - del_start) % 3:
                # Non-codon sized overlap with area of interest
                return

        for position in range(del_start, del_end):
            sequence[position - start] = ""
    elif (len(indel.what) % 3 == 0) or not options.whole_codon_indels_only:
        # parse_indel assumes that the insertion is always the first possible
        # base when multiple positions are possible. As a consequence, the
        # position may be before start, with the rest of the bases overlapping
        # the current sequence. For example:
        #  ref = ATTT
        #  alt = ATTTT
        # It is assumed that the insertion (_) happened thus:
        #  interpretation = A_TTT
        if indel.pos >= start:
            sequence[indel.pos - start] += indel.what


def filter_vcfs(genotype, contig, start, end):
    if contig in genotype.contigs:
        parser = pysam.asVCF()
        # This raises a ValueError if the VCF does not
        # contain any entries for the specified contig.
        for vcf in genotype.fetch(contig, start, end, parser=parser):
            if vcf.filter in ("PASS", "."):
                yield vcf


def build_region(options, genotype, bed):
    # Note that bed.end is a past-the-end coordinate
    start = max(0, bed.start - options.padding)

    indels = []
    sequence = ["N"] * (bed.end - start)
    for vcf in filter_vcfs(genotype, bed.contig, start, bed.end):
        if vcfwrap.is_indel(vcf):
            indels.append(vcf)
        else:
            add_snp(options, vcf, vcf.pos - start, sequence)

    if not options.ignore_indels:
        for vcf in indels:
            add_indel(options, bed, vcf, sequence)

    offset = bed.start - start
    length = bed.end - bed.start
    truncated = sequence[offset:offset + length]

    # Discard insertions after the last position
    truncated[-1] = truncated[-1][:1]

    return "".join(truncated)


def build_regions(options, genotype, beds, reverse_compl):
    for bed in beds:
        sequence = build_region(options, genotype, bed)
        if reverse_compl:
            sequence = sequences.reverse_complement(sequence)
        yield sequence


def build_genes(options, genotype, regions):
    def keyfunc(bed):
        return (bed.contig, bed.name, bed.start)
    regions.sort(key=keyfunc)

    for (gene, beds) in itertools.groupby(regions, lambda x: x.name):
        beds = split_beds(beds)
        reverse_compl = False
        if any((bed.strand == "-") for bed in beds):
            assert all((bed.strand == "-") for bed in beds)

            beds.reverse()
            reverse_compl = True

        fragments = build_regions(options, genotype, beds, reverse_compl)
        yield (gene, fragments)


def genotype_genes(options, intervals, genotype):
    for (_, beds) in sorted(intervals.items()):
        for (name, fragments) in build_genes(options, genotype, beds):
            print(">%s" % (name,))

            sequence = ""
            for fragment in fragments:
                sequence = flush_fasta(sequence + fragment)

            if sequence:
                print(sequence)

    return 0


###############################################################################
###############################################################################

def read_intervals(filename):
    with open(filename) as bed_file:
        intervals = text.parse_lines_by_contig(bed_file, BEDRecord)

        for (key, beds) in intervals.iteritems():
            bed_tuples = []
            for bed in beds:
                if len(bed) < 6:
                    sys.stderr.write(("ERROR: Invalid BED record '%r', must "
                                      "have at least 6 fields ...\n") % (bed,))
                    return None

                bed_tuples.append(bed)
            intervals[key] = bed_tuples

    return intervals


def parse_intervals(genotype):
    records = {}
    for line in genotype.header:
        match = _VCF_DICT.match(line)
        if not match:
            continue

        key, values = match.groups()
        if key == "contig":
            values = dict(pair.split("=", 1) for pair in values.split(","))
            record = BEDRecord()
            record.contig = values["ID"]
            record.start = 0
            record.end = int(values["length"])
            record.name = record.contig
            record.strand = "+"

            records[record.contig] = [record]

    if not records:
        sys.stderr.write("ERROR: List of contigs not found in VCF header; "
                         "specifying --intervals is required!\n")
        return None

    return records


def check_nth_sample(options, genotype):
    parser = pysam.asVCF()

    for contig in genotype.contigs:
        for record in text.parse_lines(genotype.fetch(contig), parser):
            if len(record) <= options.nth_sample:
                sys.stderr.write("ERROR: Sample %i selected with --nth-sample,"
                                 " but file only contains %i sample(s)!\n"
                                 % (options.nth_sample + 1, len(record)))
                return False
            return True
    return True


def main(argv):
    prog = "paleomix vcf_to_fasta"
    usage = "%s [options] --genotype in.vcf --intervals in.bed" % (prog,)

    parser = argparse.ArgumentParser(prog=prog, usage=usage)
    parser.add_argument("--genotype", required=True, metavar="VCF",
                        help="Tabix indexed VCF file; by default the first "
                             "sample is used in multi-sample VCFs. Use "
                             "--nth-sample option to select another sample.")
    parser.add_argument("--nth-sample", default=1, type=int, metavar="NTH",
                        help="Use Nth sample from the VCF, with the first "
                             "sample numbered '1' [default: %(default)s].")
    parser.add_argument("--intervals", metavar="BED",
                        help="Six column BED file; sequences on the same "
                             "contig with the same name are assumed to "
                             "represent the same gene, and are merged into a "
                             "single contiguous FASTA sequence.")
    parser.add_argument("--padding", type=int, default=10,
                        help="Number of bases to expand intervals, when "
                             "checking for adjacent indels [%(default)s]")
    parser.add_argument("--whole-codon-indels-only",
                        action="store_true", default=False,
                        help="If true, only indels where (length %% 3) == 0 "
                             "are retained [%(default)s]")
    parser.add_argument("--ignore-indels",
                        action="store_true", default=False,
                        help="Do not include indels generated FASTA "
                             "sequence [%(default)s].")

    opts = parser.parse_args(argv)

    print("Running vcf_to_fasta", end="", file=sys.stderr)
    if opts.whole_codon_indels_only:
        print(", assuming sequences represents CDS", end="", file=sys.stderr)
    print(" ...", file=sys.stderr)

    if not os.path.exists(opts.genotype):
        sys.stderr.write("ERROR: VCF file does not exist.\n")
        return 1
    elif not os.path.exists(opts.genotype + ".tbi"):
        sys.stderr.write("ERROR: VCF file not tabix indexed.\n")
        sys.stderr.write("       To index, run \"tabix -p vcf <filename>\".\n")
        return 1
    elif opts.nth_sample < 1:
        sys.stderr.write("ERROR: --nth-sample uses 1-based offsets, zero and\n")
        sys.stderr.write("       negative values are not allowed!\n")
        return 1

    # Relevant VCF functions uses zero-based offsets
    opts.nth_sample -= 1

    genotype = pysam.Tabixfile(opts.genotype)

    if opts.intervals is None:
        intervals = parse_intervals(genotype)
    else:
        intervals = read_intervals(opts.intervals)

    if intervals is None:
        return 1

    if not check_nth_sample(opts, genotype):
        return 1

    return genotype_genes(opts, intervals, genotype)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
