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

# BAM flags as defined in the BAM specification
BAM_SUPPLEMENTARY_ALIGNMENT = 0x800
BAM_PCR_DUPLICATE = 0x400
BAM_QUALITY_CONTROL_FAILED = 0x200
BAM_SECONDARY_ALIGNMENT = 0x100
BAM_IS_LAST_SEGMENT = 0x80
BAM_IS_FIRST_SEGMENT = 0x40
BAM_NEXT_IS_REVERSED = 0x20
BAM_READ_IS_REVERSED = 0x10
BAM_NEXT_IS_UNMAPPED = 0x8
BAM_READ_IS_UNMAPPED = 0x4
BAM_PROPER_SEGMENTS = 0x2
BAM_SEGMENTED = 0x1

# Default filters when processing reads
EXCLUDED_FLAGS = \
    BAM_SUPPLEMENTARY_ALIGNMENT | \
    BAM_PCR_DUPLICATE | \
    BAM_QUALITY_CONTROL_FAILED | \
    BAM_SECONDARY_ALIGNMENT | \
    BAM_READ_IS_UNMAPPED


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
    records (see e.g. paleomix.common.bedtools), with these properties:
      - contig: Name of the contig in the BED file
      - start: 0-based offset for the start of the region
      - end: 1-based offset (i.e. past-the-end) of the region
      - name: The name of the region
    """

    def __init__(self, handle, regions=None, exclude_flags=EXCLUDED_FLAGS):
        """
          - handle: BAM file handle (c.f. module 'pysam')
          - regions: List of BED-like regions (see above)
        """
        self._handle = handle
        self._regions = regions
        self._excluded = exclude_flags

    def __iter__(self):
        if self._regions:
            for region in self._regions:
                records = self._handle.fetch(region.contig,
                                             region.start,
                                             region.end)
                records = self._filter(records)

                tid = self._handle.gettid(region.contig)
                yield _BAMRegion(tid, records,
                                 region.name,
                                 region.start,
                                 region.end)
        else:
            def _by_tid(record):
                """Group by reference ID."""
                return record.tid

            # Save a copy, as these are properties generated upon every access!
            names = self._handle.references
            lengths = self._handle.lengths
            records = self._filter(self._handle)
            records = itertools.groupby(records, key=_by_tid)

            for (tid, items) in records:
                if tid >= 0:
                    name = names[tid]
                    length = lengths[tid]
                else:
                    name = length = None

                yield _BAMRegion(tid, items, name, 0, length)

    def _filter(self, records):
        """Filters records by flags, if 'exclude_flags' is set."""
        if self._excluded:
            pred = lambda record: not record.flag & self._excluded
            return itertools.ifilter(pred, records)
        return records


class _BAMRegion(object):
    """Implements iteration over sites in a BAM file. It is assumed that the
    BAM file is sorted, and that the input records are from one contig.
    """

    def __init__(self, tid, records, name, start, end):
        self._records = records
        self.tid = tid
        self.name = name
        self.start = start
        self.end = end

    def __iter__(self):
        def _by_pos(record):
            """Group by position."""
            return record.pos

        for group in itertools.groupby(self._records, _by_pos):
            yield group
