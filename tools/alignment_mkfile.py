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
#"""


_FILENAME = "SampleSheet.csv"
                     
def read_alignment_records(root):
    with open(os.path.join(root, _FILENAME)) as records:
        header = records.readline().strip().split(",")
        for line in records:
            yield dict(zip(header, line.strip().split(",")))


def main(argv):
    lines = []
    for root in argv:
        for record in read_alignment_records(root):
            record["Lane"] = int(record["Lane"])
            line = [record["SampleID"],
                    record["SampleID"],
                    record["Index"],
                    record["FCID"],
                    "Illumina",
                    os.path.join(root, "%(SampleID)s_%(Index)s_L%(Lane)03i_R{Pair}_*.fastq.gz" % record)]
            lines.append(line)
    lines.sort(key = lambda v: map(str.lower, v))
    lines.insert(0, ["Name", "Sample", "Library", "Barcode", "Platform", "Path"])

    print _HEADER % datetime.datetime.now().isoformat()
    for line in padded_table(lines):
        print line

    if not argv:
        ui.print_info("No directories specified, empty table printed:", file = sys.stderr)
        ui.print_info("\tUsage: %s [directory ...]" % sys.argv[0], file = sys.stderr)
        ui.print_info("Each directory must contain a '%s' file." % _FILENAME, file = sys.stderr)
    else:
        ui.print_info("Makefile printed. Please check for correctness and update Path column before running pipeline.", file = sys.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
