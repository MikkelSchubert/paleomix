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

from collections.abc import Generator, Iterable, Mapping
from typing import Any, TypeVar

from pysam import AlignmentFile

from paleomix.common.fileutils import open_rt
from paleomix.common.utilities import TotallyOrdered

T = TypeVar("T")


class BEDError(RuntimeError):
    pass


class BEDRecord(TotallyOrdered):
    """Class for parsing and representing a BED record.

    The class has the following properties:
       .contig -> str
       .start -> int (0-based)
       .end -> int (1-based)
       .name -> str
       .score -> int
       .strand -> '+' or '-'

    Fields past these 6 are ignored when parsing BED records.
    """

    __slots__ = ["contig", "end", "name", "score", "start", "strand"]

    def __init__(
        self,
        contig: str,
        start: int,
        end: int,
        name: str | None = None,
        score: int | None = None,
        strand: str | None = None,
    ) -> None:
        self.contig = contig
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand

        if not contig:
            raise ValueError("contig is blank")
        elif strand not in (None, "+", "-"):
            raise ValueError(f"invalid strand {strand!r}")
        elif not (0 <= start < end):
            raise ValueError("invalid start/end coordinates")

    @classmethod
    def parse(cls, line: str, _len: int | None = None) -> BEDRecord:
        """Constructs a BED record from a line of text. The length of the
        object matches the number of columns in the input line; in the case
        incompatible values, a BEDError exception is raised.
        """
        fields = line.rstrip("\r\n").split("\t")
        length = len(fields)
        if length < 3:
            raise BEDError("invalid BED record; not enough columns")
        elif not fields[0]:
            raise BEDError("contig is blank")
        elif length >= 6 and fields[5] not in "+-":
            raise BEDError("strand must be + or -")

        return BEDRecord(
            contig=fields[0],
            start=cls._parse_field(1, int, fields),
            end=cls._parse_field(2, int, fields),
            name=fields[3] if length >= 4 else None,
            score=cls._parse_field(4, int, fields) if length >= 5 else None,
            strand=fields[5] if length >= 6 else None,
        )

    @staticmethod
    def _parse_field(column: int, typ: type[T], values: list[Any]) -> T:
        try:
            return typ(values[column])
        except ValueError as error:
            raise BEDError(
                "Expected {} in column {} but found {!r}".format(
                    typ.__name__, column + 1, values[column]
                )
            ) from error

    def __str__(self) -> str:
        values = [self.contig, self.start, self.end, self.name, self.score, self.strand]
        while values and values[-1] is None:
            values.pop()

        length = len(values)
        # Default score if strand is set
        if length >= 5 and values[4] is None:
            values[4] = "0"

        # Default name if strand or score is set
        if length >= 4 and values[3] is None:
            values[3] = ""

        return "\t".join(map(str, values))

    def __repr__(self) -> str:
        keys = ("contig", "start", "end", "name", "score", "strand")
        values = [self.contig, self.start, self.end, self.name, self.score, self.strand]
        while values and values[-1] is None:
            values.pop()

        return "BEDRecord({})".format(
            ", ".join(f"{key}={value!r}" for key, value in zip(keys, values))
        )

    def __lt__(self, obj: object) -> bool:
        if not isinstance(obj, BEDRecord):
            return NotImplemented

        bed_1 = (self.contig, self.start, self.end, self.name, self.score, self.strand)
        bed_2 = (obj.contig, obj.start, obj.end, obj.name, obj.score, obj.strand)

        return bed_1 < bed_2


def read_bed_file(
    filename: str,
    contigs: Mapping[str, int] | None = None,
) -> Generator[BEDRecord, None, None]:
    """Parses a (gzip/bzip2 compressed) BED file, and yields a sequence of BED records.
    Comments and empty lines are skipped. If the number of columns in the bed record is
    less than the specified ('min_columns'), a BEDError is raised. If a dictionary of
    {contig: length} is supplied then contigs/coordinates are validated.
    """

    def _error_message(line_num: int, error: str | Exception) -> BEDError:
        return BEDError(f"{filename}:{line_num}: {error}")

    with open_rt(filename) as handle:
        for line_num, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            try:
                bed = BEDRecord.parse(line)
            except (BEDError, ValueError) as error:
                raise _error_message(line_num, error) from error

            if contigs:
                contig_len = contigs.get(bed.contig)
                if contig_len is None:
                    raise _error_message(line_num, f"unknown contig {bed.contig}")
                elif bed.end > contig_len:
                    raise _error_message(line_num, "coordinates outside contig")

            yield bed


def sort_bed_by_bamfile(bamfile: AlignmentFile, regions: list[BEDRecord]) -> None:
    """Orders a set of BED regions, such that processing matches
    (as far as possible) the layout of the BAM file. This may be
    used to ensure that extraction of regions occurs (close to)
    linearly."""
    if not regions:
        return

    references = bamfile.references
    indices = dict(zip(references, range(len(references))))
    infinite = float("inf")

    def _by_bam_layout(it: BEDRecord) -> tuple[float, str, int, int]:
        return (indices.get(it.contig, infinite), it.contig, it.start, it.end)

    regions.sort(key=_by_bam_layout)


def pad_bed_records(
    records: Iterable[BEDRecord],
    padding: int,
    max_sizes: Mapping[str, int] = {},
) -> list[BEDRecord]:
    results: list[BEDRecord] = []
    for record in records:
        start = max(0, record.start - padding)
        end = record.end + padding
        end = min(end, max_sizes.get(record.contig, end))

        if start < end:
            results.append(BEDRecord(contig=record.contig, start=start, end=end))
    return results


def merge_bed_records(records: Iterable[BEDRecord]) -> list[BEDRecord]:
    records = sorted(records)
    if not records:
        return []

    last_record = BEDRecord(
        contig=records[0].contig,
        start=records[0].start,
        end=records[0].end,
    )

    results = [last_record]
    for record in records:
        if last_record.contig != record.contig or last_record.end < record.start:
            last_record = BEDRecord(
                contig=record.contig,
                start=record.start,
                end=record.end,
            )
            results.append(last_record)
        else:
            last_record.end = record.end

    return results
