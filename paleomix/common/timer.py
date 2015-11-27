#!/usr/bin/python -3
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

import sys
import time

from paleomix.common.utilities import fragment, cumsum


_DESC  = "Processed {Records} records ({Progress}) in {Time}, est. {Remaining} left. Last {RecordsDelta} records in {TimeDelta}, now at {Contig}: {Position} ..."
_FINAL = "Processed {Records} records in {Time}. Last {RecordsDelta} records in {TimeDelta} ..."

class BAMTimer:
    def __init__(self, bamfile, desc = None, step = 1e6, out = sys.stderr):
        self._bam   = None
        self._out   = out
        self._desc  = desc
        self._step  = int(step)
        self._count = 0
        self._last_count = 0
        self._last_time  = time.time()
        self._start_time = self._last_time
        self._last_fract = -1.0

        self._bam_references = None

        self._total  = 0.0
        self._counts = []
        if bamfile and bamfile.header.get("HD", {}).get("SO", "NA") == "coordinate":
            self._bam   = bamfile
            self._bam_references = self._bam.references

            lengths = bamfile.lengths
            self._total = float(sum(lengths)) or 1.0
            self._counts.append(0)
            self._counts.extend(cumsum(lengths))


    def increment(self, count = 1, read = None):
        self._count += count
        if (self._count - self._last_count) >= self._step:
            current_time = time.time()
            self._print(current_time, read)
            self._last_time  = current_time
            self._last_count = self._count
        return self


    def finalize(self):
        self._print(time.time(), None)


    def _print(self, current_time, read):
        desc = _FINAL
        contig, position, progress, remaining = "NA", "NA", "NA", "NA"
        if read and not read.is_unmapped and self._bam:
            fraction = ((read.pos + self._counts[read.tid]) / self._total)
            if fraction >= self._last_fract:
                self._last_fract = fraction
                contig   = self._bam_references[read.tid]
                position = self._format_int(read.pos + 1)
                progress = "%.2f%%" % (fraction * 100,)

                current_running = current_time - self._start_time
                remaining = self._format_time(current_running / fraction - current_running)
                desc      = _DESC
            else:
                print("File appears to be unsorted, cannot estimate progress ...", file = self._out)
                self._bam = None

        if self._desc:
            print("%s: " % self._desc, end = "", file = self._out)

        print(desc.format(Records    = self._format_int(self._count),
                          RecordsDelta = self._format_int(self._count - self._last_count),
                          Time       = self._format_time(current_time - self._start_time),
                          TimeDelta  = self._format_time(current_time - self._last_time),
                          Contig     = contig,
                          Position   = position,
                          Progress   = progress,
                          Remaining  = remaining),
            file = self._out)


    def _format_time(self, ftime):
        utc = time.gmtime(ftime)
        return "%02i:%02i:%02is" % (utc.tm_hour, utc.tm_min, utc.tm_sec)

    def _format_int(self, value):
        return (",".join(fragment(3, str(value)[::-1])))[::-1]
