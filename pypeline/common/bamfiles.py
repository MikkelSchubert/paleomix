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
import itertools


class BAMRegionsIter(object):
    """Iterates over a BAM file, yield a separate iterator for each contig
    in the BAM or region in the list of regions if these are species, which in
    turn iterates over individual positions. This allows for the following
    pattern when parsing BAM files:

    for region in BAMRegionsIter(handle):
        # Setup per region
        for (position, records) in region:
            # Setup per position
            ...
            # Teardown per position
        # Teardown per region

    The list of regions given to the iterator is expected to be in BED-like
    records (see e.g. pypeline.common.bedtools), with these properties:
      - contig: Name of the contig in the BED file
      - start: 0-based offset for the start of the region
      - end: 1-based offset (i.e. past-the-end) of the region
      - name: The name of the region
    """

    def __init__(self, handle, regions=None):
        """
          - handle: BAM file handle (c.f. module 'pysam')
          - regions: List of BED-like regions (see above)
        """
        self._handle = handle
        self._regions = regions

    def __iter__(self):
        if self._regions:
            for region in self._regions:
                records = self._handle.fetch(region.contig,
                                             region.start,
                                             region.end)

                tid = self._handle.gettid(region.contig)
                yield _BAMRegion(tid, records,
                                 region.name,
                                 region.start,
                                 region.end)
        else:
            def _by_tid(record):
                return record.tid

            # Save a copy, as these are properties generated upon every access!
            names = self._handle.references
            lengths = self._handle.lengths
            for (tid, items) in itertools.groupby(self._handle, key=_by_tid):
                name = names[tid]
                length = lengths[tid]

                yield _BAMRegion(tid, items, name, 0, length)


class _BAMRegion(object):
    """Implements iteration over sites in a BAM file. It is assumed that the
    BAM file is sorted, and that the input records are from one contig."""

    def __init__(self, tid, records, name, start, end):
        self._records = records
        self.tid = tid
        self.name = name
        self.start = start
        self.end = end

    def __iter__(self):
        def _by_pos(record):
            return record.pos

        for group in itertools.groupby(self._records, _by_pos):
            yield group
