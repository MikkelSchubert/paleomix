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
import os
import sys
import datetime

import pypeline.ui as ui
from pypeline.common.text import padded_table

_TRIM_PIPELINE = (os.path.basename(sys.argv[0]) == "trim_pipeline")

def _print_header(timestamp, full_mkfile = True):
    print """# -*- mode: Yaml; -*-
# Timestamp: %s
#
# Default options.
# Can also be specific for a set of samples, libraries, and lanes,
# by including the "Options" hierarchy at the same level as those
# samples, libraries, or lanes below. This does not include
# "Features", which may only be specific globally.
Options:
  # Sequencing platform, see SAM/BAM reference for valid values
  Platform: Illumina
  # Quality offset for PHRED scores, either 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3+ / 1.5+)
  # For Bowtie2 it is also possible to specify 'Solexa', to handle reads on the Solexa scale.
  # This is used during adapter-trimming (AdapterRemoval) and sequence alignment (BWA/Bowtie2)
  QualityOffset: 33
  # Split a lane into multiple entries, one for each (pair of) file(s) found using the search-
  # string specified for a given lane. Each lane is named by adding a number to the end of the
  # given barcode.
  SplitLanesByFilenames: no
  # Compression format used when storing FASTQ files (either 'gz' for GZip or 'bz2' for BZip2)
  CompressionFormat: gz

  """ % timestamp

    if full_mkfile:
        print """  # Settings for aligners supported by the pipeline
  AdapterRemoval:
    # Which version of AdapterRemoval to use ('v1.4' or 'v1.5+')
    Version: 1.4

  Aligners:
    # Choice of aligner software to use, either "BWA" or "Bowtie2"
    Program: BWA

    # Settings for mappings performed using BWA
    BWA:
      # Filter hits with a mapping quality (PHRED) below this value
      MinQuality: 0
      # Should be disabled ("no") for aDNA alignments, as post-mortem localizes
      # to the seed region, which BWA expects to have few errors. Sets "-l".
      # See Schubert et al. 2012: http://pmid.us/22574660
      UseSeed:    yes
      # Additional command-line options may be specified for the "aln" call(s), as
      # described below for Bowtie2.

    # Settings for mappings performed using Bowtie2
    Bowtie2:
      # Filter hits with a mapping quality (PHRED) below this value
      MinQuality: 0
      # Examples of how to add additional command-line options
#      --trim5: 5
#      --trim3: 5
      # Note that the colon is required, even if no value is specified
      --very-sensitive:
      # Example of how to specify multiple values for an option
#     --rg:
#       - CN:SequencingCenterNameHere
#       - DS:DescriptionOfReadGroup

  # Filter PCR duplicates
  # Collapsed reads are filtered using Martin Kirchers FilterUnique,
  # while other reads are filtered using Picard MarkDuplicates.
  PCRDuplicates: yes
  # Carry out quality base re-scaling using mapDamage (*EXPERIMENTAL*)
  RescaleQualities: no

  # Exclude any type of trimmed reads from alignment/analysis
  # All reads are processed by default.
#  ExcludeReads:
#    - Single    # Single-ended reads, or PE reads where one mate was discarded
#    - Paired    # Pair-ended reads, where both reads were retained
#    - Collapsed # Overlapping pair-ended mate reads collapsed into a single read
#    - CollapsedTruncated # Like 'Collapsed', except that the reads have been
#                           truncated due to the presence of low quality bases.
#                           AdapterRemoval 1.5+ only.

  # Optional steps to perform during processing
  # To disable all features, replace with line "Features: []"
  Features:
#    - Raw BAM        # Generate BAM from the raw libraries (no indel realignment)
                     #   Location: {Destination}/{Target}.{Genome}.bam
    - Realigned BAM  # Generate indel-realigned BAM using the GATK Indel realigner
                     #   Location: {Destination}/{Target}.{Genome}.realigned.bam
    - mapDamage      # Generate mapDamage plot for each (unrealigned) library
                     #   Location: {Destination}/{Target}.{Genome}.mapDamage/{Library}/
    - Coverage       # Generate coverage information for the raw BAM (wo/ indel realignment)
                     #   Location: {Destination}/{Target}.{Genome}.coverage
    - Depths         # Generate histogram of number of sites with a given read-depth
                     #   Location: {Destination}/{Target}.{Genome}.depths
    - Summary        # Generate target summary (uses statistics from raw BAM)
                     #   Location: {Destination}/{Target}.summary


# Map of prefixes by name, each having a Path key, which specifies the location
# of the BWA/Bowtie2 index. This path should also be the filename of the
# reference FASTA sequence, such as is the case then a index is built using
# "bwa index PATH", in which case PATH would be PATH_TO_PREFIX below.
# At least ONE prefix must be specified!
#
# One or more areas of interest (for example the exome) may be specified using
# 'AreasOfInterest': Each area has a name, and points to a bedfile containing
# the relevant regions of the genome. Depths and coverage are calculated for
# these, merged by the name of the feature specified in the BED file. If no
# names are given in the BED file, the intervals are merged by contig, and
# named after the contig with a wildcard ("*") appended.
Prefixes:
  NAME_OF_PREFIX:
    Path: PATH_TO_PREFIX
#    Label: # "mito" or "nucl"
#    AreasOfInterest:
#      NAME: PATH_TO_BEDFILE

# Prefixes can also be specified using wildcards, by adding a wildcard to the
# end of the prefix name, as shown below. The name itself is ignored, and each
# prefix is named according to the basename of the path (ie. the filename with
# the extensions removed).
#  NAME_OR_DESCRIPTION*:
#    Path: PATH_TO_PREFIXES/*.fasta
#    Label: # "mito" or "nucl"

"""


_FILENAME = "SampleSheet.csv"

def read_alignment_records(filename):

    with open(filename) as records:
        header = records.readline().strip().split(",")
        for line in records:
            yield dict(zip(header, line.strip().split(",")))


def main(argv):
    records = {}
    for root in argv:
        if os.path.isdir(root):
            filename = os.path.join(root, _FILENAME)
        else:
            root, filename = os.path.split(root)[0], root

        for record in read_alignment_records(filename):
            libraries = records.setdefault(record["SampleID"], {})
            barcodes  = libraries.setdefault(record["Index"], [])

            record["Lane"] = int(record["Lane"])
            record["Path"] = os.path.join(root, "%(SampleID)s_%(Index)s_L%(Lane)03i_R{Pair}_*.fastq.gz" % record)
            barcodes.append(record)


    _print_header(timestamp   = datetime.datetime.now().isoformat(),
                  full_mkfile = (os.path.basename(sys.argv[0]) != "trim_pipeline"))
    for (sample, libraries) in records.iteritems():
        print "%s:" % sample
        print "  %s:" % sample
        for (library, barcodes) in libraries.iteritems():
            print "    %s:" % library
            for record in barcodes:
                print "      {FCID}_{Lane}: {Path}".format(**record)
            print
        print

    if not argv:
        ui.print_info("No directories specified, empty table printed:", file = sys.stderr)
        ui.print_info("\tUsage: %s [directory ...]" % sys.argv[0], file = sys.stderr)
        ui.print_info("Each directory must contain a '%s' file." % _FILENAME, file = sys.stderr)
    else:
        ui.print_info("Makefile printed. Please check for correctness and update Path column before running pipeline.", file = sys.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
