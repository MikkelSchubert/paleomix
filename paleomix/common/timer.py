#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import sys
import time
from typing import IO, Iterator

from pysam import AlignedSegment, AlignmentFile


class BAMTimer:
    def __init__(
        self,
        handle: AlignmentFile,
        step: int = 1_000_000,
        out: IO[str] = sys.stderr,
    ):
        self._handle = handle
        self._out = out
        self._step = step
        self._count = 0
        self._last_count = 0
        self._last_time = time.monotonic()
        self._start_time = self._last_time

    def increment(self, count: int = 1) -> None:
        self._count += count
        if self._count - self._last_count >= self._step:
            current_time = time.monotonic()
            self._print(current_time)
            self._last_time = current_time
            self._last_count = self._count

    def finalize(self) -> None:
        self._print(time.monotonic())

    def __iter__(self) -> Iterator[AlignedSegment]:
        for record in self._handle:
            yield record
            self.increment()

        self.finalize()

    def _print(self, current_time: float) -> None:
        print(
            "Processed {:,} records in {}. Last {:,} records in {}".format(
                self._count,
                self._format_time(current_time - self._start_time),
                self._count - self._last_count,
                self._format_time(current_time - self._last_time),
            ),
            file=self._out,
        )

    def _format_time(self, secs: float) -> str:
        return time.strftime("%H:%M:%Ss", time.gmtime(secs))
