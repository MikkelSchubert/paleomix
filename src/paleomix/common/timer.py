# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import sys
import time
from collections.abc import Iterator
from typing import IO

from pysam import AlignedSegment, AlignmentFile


class BAMTimer:
    def __init__(
        self,
        handle: AlignmentFile,
        step: int = 1_000_000,
        out: IO[str] = sys.stderr,
    ) -> None:
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
