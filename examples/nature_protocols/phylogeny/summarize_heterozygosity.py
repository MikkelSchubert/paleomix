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
import sys
import pysam

from paleomix.common.vcfwrap import \
     get_ml_genotype

import paleomix.common.timer as timer


def read_bed_records(filename):
    """Reads a bed-file (i.e. for a set of regions of interest), and returns
    a sorted list containing each line as a tuple containing the contig name,
    the start position, and the end position."""
    regions = []
    bed_parser = pysam.asBed()
    with open(filename) as bed_file:
        for line in bed_file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            regions.append(bed_parser(line, len(line)))
    return regions


def select_vcf_records(bed_records, vcf_records):
    """Returns an iterable of VCF records, corresponding to the contents of each
    region specified by the BED records. Records are returned at most once, even
    if covered by multiple BED records."""
    contigs    = frozenset(vcf_records.contigs)
    vcf_parser = pysam.asVCF()

    # Timer class used processing progress; meant primarily for BAM files
    progress   = timer.BAMTimer(None)

    # Cache of positions observed for this contig, to prevent returning
    # positions in overlapping regions multiple times
    contig_cache = None
    contig_cache_name = None

    for bed in sorted(bed_records):
        if bed.contig not in contigs:
            # Skip contigs for which no calls have been made (e.g. due to
            # low coverage. Otherwise Pysam raises an exception.
            continue
        elif contig_cache_name != bed.contig:
            # Reset cache per contig, to save memory
            contig_cache = set()
            contig_cache_name = bed.contig

        for record in vcf_records.fetch(bed.contig, bed.start, bed.end, parser = vcf_parser):
            progress.increment()

            if record.pos in contig_cache:
                # We've already reported this VCF record
                continue

            contig_cache.add(record.pos)
            # Skip records filtered by VCF_filter
            if record.filter in ('.', "PASS"):
                yield record
    progress.finalize()


def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <BED-file> <VCF.bgz>\n")
        return 1

    sites                 = 0
    sites_non_ref         = 0
    sites_homo_non_ref    = 0
    sites_het_one_non_ref = 0
    sites_het_two_non_ref = 0

    vcf_records = pysam.Tabixfile(argv[1])
    bed_records = read_bed_records(argv[0])

    for record in select_vcf_records(bed_records, vcf_records):
        if record.alt != '.':
            # Get the most likely diploid genotype
            nt_a, nt_b = get_ml_genotype(record)
            if (nt_a, nt_b) == ('N', 'N'):
                # Skip sites with no most likely genotype
                continue

            sites += 1
            sites_non_ref += 1
            if nt_a == nt_b:
                sites_homo_non_ref += 1
            elif record.ref not in (nt_a, nt_b):
                sites_het_two_non_ref += 1
            else:
                sites_het_one_non_ref += 1
        else:
            # Heterozygous for the reference allele
            sites += 1

    print
    print "%i sites kept after filtering:" % (sites,)
    print " % 10i homozygous sites containing the reference allele (%.2f%%)" % (sites - sites_non_ref, 100.0 * (sites - sites_non_ref) / float(sites))
    print " % 10i heterozygous sites containing the reference and a non-reference allele (%.2f%%)" % (sites_het_one_non_ref, (100.0 * sites_het_one_non_ref) / sites)
    print " % 10i homozygous sites containing a single non-reference allele (%.2f%%)" % (sites_homo_non_ref, (100.0 * sites_homo_non_ref) / sites)
    print " % 10i heterozygous sites containing two different non-reference alleles (%.2f%%)" % (sites_het_two_non_ref, (100.0 * sites_het_two_non_ref) / sites)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
