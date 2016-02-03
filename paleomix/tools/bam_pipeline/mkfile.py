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

import os
import sys
import glob
import datetime
from optparse import OptionParser

from paleomix.common.console import \
    print_info, \
    print_err


_TEMPLATE_TOP = \
    """# -*- mode: Yaml; -*-
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
  # Quality offset for Phred scores, either 33 (Sanger/Illumina 1.8+)
  # or 64 (Illumina 1.3+ / 1.5+). For Bowtie2 it is also possible to
  # specify 'Solexa', to handle reads on the Solexa scale. This is
  # used during adapter-trimming and sequence alignment
  QualityOffset: 33
  # Split a lane into multiple entries, one for each (pair of) file(s)
  # found using the search-string specified for a given lane. Each
  # lane is named by adding a number to the end of the given barcode.
  SplitLanesByFilenames: yes
  # Compression format for FASTQ reads; 'gz' for GZip, 'bz2' for BZip2
  CompressionFormat: bz2

  # Settings for trimming of reads, see AdapterRemoval man-page
  AdapterRemoval:
     # Adapter sequences, set and uncomment to override defaults
#     --adapter1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#     --adapter2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
     # Some BAM pipeline defaults differ from AR defaults;
     # To override, change these value(s):
     --mm: 3
     --minlength: 25
     # Extra features enabled by default; change 'yes' to 'no' to disable
     --collapse: yes
     --trimns: yes
     --trimqualities: yes
"""

_TEMPLATE_BAM_OPTIONS = \
    """  # Settings for aligners supported by the pipeline
  Aligners:
    # Choice of aligner software to use, either "BWA" or "Bowtie2"
    Program: BWA

    # Settings for mappings performed using BWA
    BWA:
      # One of "backtrack", "bwasw", or "mem"; see the BWA documentation
      # for a description of each algorithm (defaults to 'backtrack')
      Algorithm: backtrack
      # Filter aligned reads with a mapping quality (Phred) below this value
      MinQuality: 0
      # Filter reads that did not map to the reference sequence
      FilterUnmappedReads: yes
      # Should be disabled ("no") for aDNA alignments, as post-mortem
      # localizes to the seed region, which BWA expects to have few
      # errors (sets "-l"). See http://pmid.us/22574660
      UseSeed:    yes
      # Additional command-line options may be specified for the "aln"
      # call(s), as described below for Bowtie2 below.

    # Settings for mappings performed using Bowtie2
    Bowtie2:
      # Filter aligned reads with a mapping quality (Phred) below this value
      MinQuality: 0
      # Filter reads that did not map to the reference sequence
      FilterUnmappedReads: yes
      # Examples of how to add additional command-line options
#      --trim5: 5
#      --trim3: 5
      # Note that the colon is required, even if no value is specified
      --very-sensitive:
      # Example of how to specify multiple values for an option
#      --rg:
#        - CN:SequencingCenterNameHere
#        - DS:DescriptionOfReadGroup

  # Mark / filter PCR duplicates. If set to 'filter', PCR duplicates are
  # removed from the output files; if set to 'mark', PCR duplicates are
  # flagged with bit 0x400, and not removed from the output files; if set to
  # 'no', the reads are assumed to not have been amplified. Collapsed reads
  # are filtered using the command 'paleomix rmdup_duplicates', while "normal"
  # reads are filtered using Picard MarkDuplicates.
  PCRDuplicates: filter

  # Carry out quality base re-scaling of libraries using mapDamage
  # This will be done using the options set for mapDamage below
  RescaleQualities: no

  # Command-line options for mapDamage; note that the long-form
  # options are expected; --length, not -l, etc. Uncomment the
  # "mapDamage" line adding command-line options below.
  mapDamage:
    # By default, the pipeline will downsample the input to 100k hits
    # when running mapDamage; remove to use all hits
    --downsample: 100000

  # Set to 'yes' exclude a type of trimmed reads from alignment / analysis;
  # possible read-types reflect the output of AdapterRemoval
  ExcludeReads:
    Single: no              # Single-ended reads / Orphaned paired-ended reads
    Paired: no              # Paired ended reads
    Singleton: no           # Paired reads for which the mate was discarded
    Collapsed: no           # Overlapping paired-ended reads collapsed into a
                            # single sequence by AdapterRemoval
    CollapsedTruncated: no  # Like 'Collapsed', except that the reads
                            # truncated due to the presence ambigious
                            # bases or low quality bases at read termini.

  # Optional steps to perform during processing
  Features:
    RawBAM: no          # Generate BAM from the raw libraries (no indel realignment)
                        #   Location: {Destination}/{Target}.{Genome}.bam
    RealignedBAM: yes   # Generate indel-realigned BAM using the GATK Indel realigner
                        #   Location: {Destination}/{Target}.{Genome}.realigned.bam
    mapDamage: yes      # Generate mapDamage plot for each (unrealigned) library
                        #   Location: {Destination}/{Target}.{Genome}.mapDamage/{Library}/
    Coverage: yes       # Generate coverage information for the raw BAM (wo/ indel realignment)
                        #   Location: {Destination}/{Target}.{Genome}.coverage
    Depths: yes         # Generate histogram of number of sites with a given read-depth
                        #   Location: {Destination}/{Target}.{Genome}.depths
    Summary: yes        # Generate summary table for each target
                        #   Location: {Destination}/{Target}.summary
    DuplicateHist: no   # Generate histogram of PCR duplicates, for use with PreSeq
                        #   Location: {Destination}/{Target}.{Genome}.duphist/{Library}/
"""

_TEMPLATE_PREFIXES = """
# Map of prefixes by name, each having a Path key, which specifies the
# location of the BWA/Bowtie2 index, and optional label, and an option
# set of regions for which additional statistics are produced.
Prefixes:
  # Name of the prefix; is used as part of the output filenames
  NAME_OF_PREFIX:
    # Path to .fasta file containg a set of reference sequences.
    Path: PATH_TO_PREFIX

    # Label for prefix: One of nuclear, mitochondrial, chloroplast,
    # plasmid, bacterial, or viral. Is used in the .summary files.
#    Label: ...

    # Produce additional coverage / depth statistics for a set of
    # regions defined in a BED file; if no names are specified for the
    # BED records, results are named after the chromosome / contig.
#    RegionsOfInterest:
#      NAME: PATH_TO_BEDFILE
"""

_TEMPLATE_SAMPLES = """
# Targets are specified using the following structure:
#NAME_OF_TARGET:
#  NAME_OF_SAMPLE:
#    NAME_OF_LIBRARY:
#      NAME_OF_LANE: PATH_WITH_WILDCARDS
"""

_FILENAME = "SampleSheet.csv"


def build_makefile(add_full_options=True,
                   add_prefix_tmpl=True,
                   add_sample_tmpl=True):
    timestamp = datetime.datetime.now().isoformat()
    template_parts = [_TEMPLATE_TOP % (timestamp,)]

    if add_full_options:
        template_parts.append(_TEMPLATE_BAM_OPTIONS)

    if add_prefix_tmpl:
        template_parts.append(_TEMPLATE_PREFIXES)

    if add_sample_tmpl:
        template_parts.append(_TEMPLATE_SAMPLES)

    return "\n".join(template_parts)


def strip_comments(text):
    lines = text.split("\n")

    # Always include minimal header
    minimal_template = lines[:3]
    for line in lines[3:]:
        if not line.lstrip().startswith("#"):
            line = line.split("#", 1)[0].rstrip()

            # Avoid too many empty lines
            if line.strip() or minimal_template[-1].strip():
                minimal_template.append(line)

    return "\n".join(minimal_template)


def read_alignment_records(filename):
    with open(filename) as records:
        header = records.readline().strip().split(",")
        for line in records:
            yield dict(zip(header, line.strip().split(",")))


def parse_args(argv):
    parser = OptionParser("Usage: %prog [/path/to/SampleSheet.csv, ...]")
    parser.add_option("--minimal", default=False, action="store_true",
                      help="Strip comments from makefile template.")

    return parser.parse_args(argv)


def select_path(path):
    has_r1 = bool(glob.glob(path.format(Pair=1)))
    has_r2 = bool(glob.glob(path.format(Pair=2)))

    if has_r1 and not has_r2:
        # Single-ended reads
        return path.format(Pair=1)
    return path


def main(argv, pipeline="bam"):
    assert pipeline in ("bam", "trim"), pipeline

    options, paths = parse_args(argv)
    records = {}
    for root in paths:
        if os.path.isdir(root):
            filename = os.path.join(root, _FILENAME)
        else:
            root, filename = os.path.split(root)[0], root

        if not os.path.exists(filename):
            print_err("ERROR: Could not find SampleSheet file: %r" % filename)
            return 1

        for record in read_alignment_records(filename):
            libraries = records.setdefault(record["SampleID"], {})
            barcodes = libraries.setdefault(record["Index"], [])

            record["Lane"] = int(record["Lane"])
            path = "%(SampleID)s_%(Index)s_L%(Lane)03i_R{Pair}_*.fastq.gz" \
                % record
            record["Path"] = select_path(os.path.join(root, path))
            barcodes.append(record)

    template = build_makefile(add_full_options=(pipeline == "bam"))
    if options.minimal:
        template = strip_comments(template)

    print(template)

    for (sample, libraries) in records.iteritems():
        print("%s:" % sample)
        print("  %s:" % sample)
        for (library, barcodes) in libraries.iteritems():
            print("    %s:" % library)
            for record in barcodes:
                print("      {FCID}_{Lane}: {Path}".format(**record))
            print()
        print()

    if argv:
        print_info("Automatically generated makefile printed.\n"
                   "Please check for correctness before running pipeline.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
