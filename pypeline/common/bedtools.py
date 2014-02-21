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
import pysam

import pypeline.common.text as text
import pypeline.common.fileutils as fileutils


def read_bed_file(filename):
    """Parses a (gzip/bzip2 compressed) BED file, and yields
    a sequence of records. Comments and empty lines are skipped."""
    handle = None
    try:
        handle = fileutils.open_ro(filename)
        parser = pysam.asBed()

        for record in text.parse_lines(handle, parser):
            # Force evaluation of (lazily parsed) properties
            _ = record.start
            _ = record.end

            yield record

    finally:
        if handle:
            handle.close()


def sort_bed_by_bamfile(bamfile, regions):
    """Orders a set of BED regions, such that processing matches
    (as far as possible) the layout of the BAM file. This may be
    used to ensure that extraction of regions occurs (close to)
    linearly."""
    if not regions:
        return

    indices = dict(zip(bamfile.references,
                   xrange(len(bamfile.references))))

    def _by_bam_layout(region):
        return (indices[region.contig], region.start, region.end)
    regions.sort(key=_by_bam_layout)
