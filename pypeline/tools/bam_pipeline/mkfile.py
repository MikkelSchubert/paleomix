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


_HEADER = \
"""# Timestamp: %s
# 
#
#
# Columns:
#     Name:     Name of the target. Final BAM filenames will consist of this
#               name and the name of the BWA prefix joined by a dot. For
#               example, given the name 'foobar' and the BWA prefix 
#               'path/to/indices/hg19', the final BAM files will be named
#               'foobar.hg19.bam' and 'foobar.hg19.unaligned.bam'.
#     Sample:
#     Library:  
#     Barcode:  The Sample, Library, and Barcode (or other unique ID) of the
#               lane. This is saved in the 'SM', 'LB' and 'PU'. 
#     Platform: The sequencing platform used to generate this data. 
#     Path:     Path to FASTQ files. This uses the Python 'glob' module to find
#               files, and may therefore contain wildcards (*, [...]). For PE
#               reads, replace the mate number with "{Pair}". 
#               For example
#                  /path/to/files/reads_R{Pair}_*.fastq.gz
#
# Data is organized by Name/Sample/Library/Barcode. Therefore, this combination
# of columns MUST be unique. To prevent duplicate removal being done on multiple 
# different samples, it is enforced that all libraries are assosiated with ONE 
# sample only for any one target. Finally, any input file may only be specified
# ONCE. The pipeline will not run if any of these assumptions are violated.
#
#
"""

_HEADER = \
"""# Timestamp: %s
#
# Default options
Options:
  # Sequencing platform, see SAM/BAM reference for valid values
  Platform: Illumina
  
  # Use seed during sequence alignment
  BWA_UseSeed: yes
  # Filter hits with a mappign quality (PHRED) below this value
  BWA_MinQuality: 0

  # Filter PCR duplicates
  PCRDuplicates: yes

  # Exclude any type of trimmed reads from alignment/analysis
  # All reads are processed by default.
#  ExcludeReads:
#  - Single    # Single-ended reads, or PE reads where one mate was discarded
#  - Paired    # Pair-ended reads, where both reads were retained
#  - Collapsed # Overlapping pair-ended mate reads collapsed into a single read  # ExcludeReads: 


# Prefixes:
#    - NAME_OF_PREFIX: PATH_TO_PREFIX
#      Label: # mito or nucl
#
#
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


    print _HEADER % datetime.datetime.now().isoformat()
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
