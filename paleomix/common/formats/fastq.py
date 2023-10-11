#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import sys
from enum import Enum
from typing import IO, TYPE_CHECKING, Iterable, Iterator

from paleomix.common.fileutils import open_rt
from paleomix.common.formats._common import FormatError
from paleomix.common.utilities import Immutable, TotallyOrdered

if TYPE_CHECKING:
    from pathlib import Path


class FASTQError(FormatError):
    pass


class FASTQ(TotallyOrdered, Immutable):
    __slots__ = ["name", "meta", "sequence", "qualities"]

    name: str
    meta: str | None
    sequence: str
    qualities: str

    def __init__(
        self,
        name: str,
        meta: str | None,
        sequence: str,
        qualities: str,
    ) -> None:
        if not (name and isinstance(name, str)):
            raise FASTQError("FASTQ name must be a non-empty string")
        elif not (isinstance(meta, str) or (meta is None)):
            raise FASTQError("FASTQ meta must be a string, or None")
        elif not isinstance(sequence, str):
            raise FASTQError("FASTQ sequence must be a string")
        elif not isinstance(qualities, str):
            raise FASTQError("FASTQ qualities must be a string")
        elif len(sequence) != len(qualities):
            raise FASTQError("FASTQ sequence and qualities length differ")

        Immutable.__init__(
            self, name=name, meta=meta or "", sequence=sequence, qualities=qualities
        )

    def write(self, fileobj: IO[str] = sys.stdout) -> None:
        """Writes a FASTQ record to fileobj."""
        name = self.name
        if self.meta:
            name = f"{name} {self.meta}"

        fileobj.write(f"@{name}\n{self.sequence}\n+\n{self.qualities}\n")

    @classmethod
    def from_lines(cls, lines: Iterable[str]) -> Iterator[FASTQ]:
        """Parses FASTQ sequences found in a sequence of lines, and returns
        a tuple for each FASTQ record: ((name, meta-information), sequence)
        No assumptions are made about the line-lengths."""
        lines_iter = iter(lines)
        while True:
            try:
                header = next(lines_iter).rstrip()
                if not header:
                    break
            except StopIteration:
                break  # No sequences left

            if not header.startswith("@"):
                raise FASTQError(f"Invalid FASTQ header: {header!r}]")

            header_fields = header.split(None, 1)
            name = header_fields[0]
            meta = header_fields[1] if len(header_fields) == 2 else None

            try:
                sequence = next(lines_iter).rstrip()
                separator = next(lines_iter).rstrip()
                qualities = next(lines_iter).rstrip()
            except StopIteration as error:
                raise FASTQError(f"Partial FASTQ record: {name}") from error

            if not separator.startswith("+"):
                raise FASTQError(
                    f"Invalid FASTQ separator for {name!r}; "
                    f"expected '+', found {separator!r}"
                )
            elif len(sequence) != len(qualities):
                raise FASTQError(
                    f"Sequence length does not match qualities length for {name!r}"
                )

            yield FASTQ(
                name=name[1:], meta=meta, sequence=sequence, qualities=qualities
            )

    @classmethod
    def from_file(cls, filename: str | Path) -> Iterator[FASTQ]:
        """Reads an unindexed FASTQ file, returning a sequence of
        tuples containing the name and sequence of each entry in
        the file. The FASTQ file may be GZIP/BZ2 compressed."""
        with open_rt(filename) as handle:
            yield from FASTQ.from_lines(handle)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, FASTQ):
            return NotImplemented

        return (self.name, self.meta, self.sequence, self.qualities) < (
            other.name,
            other.meta,
            other.sequence,
            other.qualities,
        )

    def __hash__(self) -> int:
        return hash((self.name, self.meta, self.sequence, self.qualities))

    def __repr__(self) -> str:
        return "FASTQ({!r}, {!r}, {!r}, {!r})".format(
            self.name,
            self.meta,
            self.sequence,
            self.qualities,
        )


class FASTQOffsets(Enum):
    # Quality score offsets for Phred (or similar) scores in FASTQ reads (33 or 64)
    OFFSET_33 = 33
    OFFSET_64 = 64
    # Quality score found in both ranges that are unique to each offset,
    # suggesting that the list contains mixed quality offsets, or invalid data.
    BOTH = "BOTH"
    # No quality scores in the expected range of ASCII values (33 .. 105)
    MISSING = "MISSING"
    # Quality scores are in the ASCII range 59 .. 74, which could signify
    # low-quality reads with offset 64, or high-quality reads with offset 33
    AMBIGIOUS = "AMBIGIOUS"


class FASTQualities:
    """Given a set of FASTQ records, this class attempts to identify the
    offset used to encode the quality scores pf those records.

    The following constants may be returned:
    - OFFSET_33: Offset identified as being 33
    - OFFSET_64: Offset identified as being 64
    - BOTH: Both offset 33 and 64 found, mixed file? (error)
    - MISSING: No quality scores found, wrong file? (error)
    - AMBIGIOUS: Qualities could be either offset. (warning)
    """

    def __init__(self) -> None:
        self._qualities: set[str] = set()

    def update(self, record: FASTQ) -> None:
        self._qualities.update(record.qualities)

    def offsets(self) -> FASTQOffsets:
        qualities = [False] * 256
        for quality in self._qualities:
            qualities[ord(quality)] = True

        # The range of scores that can unambigiously be identified
        # as belonging to Phred scores with offset 33 or 64. Scores
        # in between could potentially signify either offset
        # See e.g. http://en.wikipedia.org/wiki/FASTQ_format#Encoding
        has_offset_33_scores = any(qualities[33:59])
        has_ambigious_scores = any(qualities[59:75])
        has_offset_64_scores = any(qualities[75:105])

        if has_offset_33_scores:
            if has_offset_64_scores:
                return FASTQOffsets.BOTH
            return FASTQOffsets.OFFSET_33
        elif has_offset_64_scores:
            return FASTQOffsets.OFFSET_64
        elif has_ambigious_scores:
            return FASTQOffsets.AMBIGIOUS

        return FASTQOffsets.MISSING
