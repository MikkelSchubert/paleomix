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
"""Tool for building histogram of PCR duplicates to be used with 'preseq'.

This allows estimation of the library complexity, and potential further gains
from sequencing the library. Unlike the tools included in 'preseq', this tool
handles collapsed reads.

Preseq is located at http://smithlabresearch.org/software/preseq/
"""
import sys
import random
import argparse
import collections

import pysam

import paleomix.common.bamfiles as bamfiles


def get_template_length(record):
    """Returns the template length of the given record for paired or collapsed
    reads; for single-ended reads, None is returned."""
    if record.is_paired:
        return record.tlen
    elif record.qname.startswith("M_"):
        return record.alen
    return None


def process_records(records, counts):
    """Processes a set of records aligned to the same positon; dulpicates are
    inferred based on the template lengths of the records, using the sequence
    length for collapsed reads (where qname.startswith("M_")). Single-ended
    reads are assumed to represent a random sampling of reads for which the
    insert size is known.
    """
    alignments = collections.defaultdict(int)
    for record in records:
        if record.is_paired:
            if (not record.is_proper_pair) or record.is_read2:
                continue

        alignment = get_template_length(record)
        alignments[alignment] += 1

    if (None in alignments) and len(alignments) > 1:
        ambigious_count = alignments.pop(None)

        # PE reads are assummed to represent a random sample of PE / collapsed
        # reads, and distributed randomly across these to approximate this
        keys = tuple(alignments)
        for _ in xrange(ambigious_count):
            key = random.choice(keys)
            alignments[key] += 1

    for count in alignments.itervalues():
        counts[count] += 1


def parse_args(argv):
    prog = "paleomix duphist"
    usage = "%s sorted.bam > out.histogram" % (prog,)
    parser = argparse.ArgumentParser(prog=prog, usage=usage)
    parser.add_argument("bamfile", help="Sorted BAM file.")

    return parser.parse_args(argv)


def main(argv):
    """Main function; takes a list of arguments equivalent to sys.argv[1:]."""
    args = parse_args(argv)

    # Default filters, excepting that PCR duplicates are not filtered
    mask = bamfiles.EXCLUDED_FLAGS & ~bamfiles.BAM_PCR_DUPLICATE
    counts = collections.defaultdict(int)
    with pysam.Samfile(args.bamfile) as handle:
        for region in bamfiles.BAMRegionsIter(handle, exclude_flags=mask):
            for (_, records) in region:
                process_records(records, counts)

    for (key, count) in sorted(counts.iteritems()):
        print "%i\t%i" % (key, count)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
