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
      # May be disabled ("no") for aDNA alignments, as post-mortem damage
      # localizes to the seed region, which BWA expects to have few
      # errors (sets "-l"). See http://pmid.us/22574660
      UseSeed: yes
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
    # Exclude single-end reads (yes / no)?
    Single: no
    # Exclude non-collapsed paired-end reads (yes / no)?
    Paired: no
    # Exclude paired-end reads for which the mate was discarded (yes / no)?
    Singleton: no
    # Exclude overlapping paired-ended reads collapsed into a single sequence
    # by AdapterRemoval (yes / no)?
    Collapsed: no
    # Like 'Collapsed', but only for collapsed reads truncated due to the
    # presence of ambiguous or low quality bases at read termini (yes / no).
    CollapsedTruncated: no

  # Optional steps to perform during processing.
  Features:
    # Generate BAM without realignment around indels (yes / no)
    RawBAM: no
    # Generate indel-realigned BAM using the GATK Indel realigner (yes / no)
    RealignedBAM: yes
    # To disable mapDamage, write 'no'; to generate basic mapDamage plots,
    # write 'plot'; to build post-mortem damage models, write 'model',
    # and to produce rescaled BAMs, write 'rescale'. The 'model' option
    # includes the 'plot' output, and the 'rescale' option includes both
    # 'plot' and 'model' results. All analyses are carried out per library.
    mapDamage: plot
    # Generate coverage information for the raw BAM (wo/ indel realignment).
    # If one or more 'RegionsOfInterest' have been specified for a prefix,
    # additional coverage files are generated for each alignment (yes / no)
    Coverage: yes
    # Generate histogram of number of sites with a given read-depth, from 0
    # to 200. If one or more 'RegionsOfInterest' have been specified for a
    # prefix, additional histograms are generated for each alignment (yes / no)
    Depths: yes
    # Generate summary table for each target (yes / no)
    Summary: yes
    # Generate histogram of PCR duplicates, for use with PreSeq (yes / no)
    DuplicateHist: no
"""

_TEMPLATE_PREFIXES = """
# Map of prefixes by name, each having a Path key, which specifies the
# location of the BWA/Bowtie2 index, and optional label, and an option
# set of regions for which additional statistics are produced.
Prefixes:
  # Replace 'NAME_OF_PREFIX' with name of the prefix; this name
  # is used in summary statistics and as part of output filenames.
  NAME_OF_PREFIX:
    # Replace 'PATH_TO_PREFIX' with the path to .fasta file containing the
    # references against which reads are to be mapped. Using the same name
    # as filename is strongly recommended (e.g. /path/to/Human_g1k_v37.fasta
    # should be named 'Human_g1k_v37').
    Path: PATH_TO_PREFIX

    # (Optional) Uncomment and replace 'PATH_TO_BEDFILE' with the path to a
    # .bed file listing extra regions for which coverage / depth statistics
    # should be calculated; if no names are specified for the BED records,
    # results are named after the chromosome / contig. Change 'NAME' to the
    # name to be used in summary statistics and output filenames.
#    RegionsOfInterest:
#      NAME: PATH_TO_BEDFILE
"""

_TEMPLATE_SAMPLES = """
# Mapping targets are specified using the following structure. Uncomment and
# replace 'NAME_OF_TARGET' with the desired prefix for filenames.
#NAME_OF_TARGET:
   #  Uncomment and replace 'NAME_OF_SAMPLE' with the name of this sample.
#  NAME_OF_SAMPLE:
     #  Uncomment and replace 'NAME_OF_LIBRARY' with the name of this sample.
#    NAME_OF_LIBRARY:
       # Uncomment and replace 'NAME_OF_LANE' with the name of this lane,
       # and replace 'PATH_WITH_WILDCARDS' with the path to the FASTQ files
       # to be trimmed and mapped for this lane (may include wildcards).
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
    results = []
    with open(filename) as records:
        line = records.readline()
        if not line:
            print_err("ERROR: Empty SampleSheet.csv file: %r"
                      % (filename,))
            return None

        header = line.strip().split(",")
        missing = set(("SampleID", "Index", "Lane", "FCID")) - set(header)
        if missing:
            print_err("ERROR: Required columns missing from SampleSheet file "
                      "%r: %s" % (filename, ", ".join(map(repr, missing))))
            return None

        for idx, line in enumerate(records, start=2):
            line = line.strip()
            if not line:
                continue

            fields = line.split(",")
            if len(fields) != len(header):
                print_err("Line %i in SampleSheet file %r does not contain "
                          "the expected number of columns; expected %i, but "
                          "found %i."
                          % (idx, filename, len(header), len(fields)))
                return None

            results.append(dict(zip(header, fields)))

    return results


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


def read_sample_sheets(filenames):
    records = {}
    for root in filenames:
        if os.path.isdir(root):
            filename = os.path.join(root, _FILENAME)
        else:
            root, filename = os.path.split(root)[0], root

        if not os.path.exists(filename):
            print_err("ERROR: Could not find SampleSheet file: %r" % filename)
            return None

        sample_sheet = read_alignment_records(filename)
        if sample_sheet is None:
            return None

        for record in sample_sheet:
            record["Lane"] = int(record["Lane"])
            path = "%(SampleID)s_%(Index)s_L%(Lane)03i_R{Pair}_*.fastq.gz" \
                % record
            record["Path"] = select_path(os.path.join(root, path))
            key = "%(FCID)s_%(Lane)s" % record

            libraries = records.setdefault(record["SampleID"], {})
            barcodes = libraries.setdefault(record["Index"], {})
            barcodes.setdefault(key, []).append(path)

    # Clean up names; generate unique names for duplicate lanes
    for libraries in records.itervalues():
        for barcodes in libraries.itervalues():
            for key, paths in barcodes.items():
                if len(paths) == 1:
                    barcodes[key] = paths[0]
                    continue

                counter = 1
                for path in paths:
                    new_key = "%s_%i" % (key, counter)

                    while new_key in barcodes:
                        counter += 1
                        new_key = "%s_%i" % (key, counter)

                    barcodes[new_key] = path

                barcodes.pop(key)

    return records


def print_samples(records):
    print()
    for (sample, libraries) in sorted(records.iteritems()):
        print("%s:" % sample)
        print("  %s:" % sample)
        for (library, barcodes) in sorted(libraries.iteritems()):
            print("    %s:" % library)
            for key, path in sorted(barcodes.iteritems()):
                print("      %s: %s" % (key, path))
            print()
        print()


def main(argv, pipeline="bam"):
    assert pipeline in ("bam", "trim"), pipeline

    options, filenames = parse_args(argv)
    records = read_sample_sheets(filenames)
    if records is None:
        return 1

    template = build_makefile(add_full_options=(pipeline == "bam"),
                              add_prefix_tmpl=(pipeline == "bam"),
                              add_sample_tmpl=not records)
    if options.minimal:
        template = strip_comments(template)

    print(template)

    print_samples(records)

    if argv:
        print_info("Automatically generated makefile printed.\n"
                   "Please check for correctness before running pipeline.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
