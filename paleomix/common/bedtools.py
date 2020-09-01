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
import copy

from paleomix.common.fileutils import open_ro
from paleomix.common.utilities import TotallyOrdered


def _strand_type(value):
    if value not in ("+", "-"):
        raise ValueError("Strand must be '+' or '-', not %r" % (value,))
    return value


_BED_DEFAULTS = ("", 0, 0, "", 0, "+")
_BED_KEYS = ("contig", "start", "end", "name", "score", "strand")
_BED_TYPES = (str, int, int, str, int, _strand_type)


class BEDError(RuntimeError):
    pass


class BEDRecord(TotallyOrdered):
    """Class for parsing and representing a BED records.

    The class has the following properties:
       .contig -> str
       .start -> int (0-based)
       .end -> int (1-based)
       .name -> str
       .score -> int
       .strand -> '+' or '-'

    These fields can also be accessed using the square brackets notation, which
    also gives access to any additional values, after the strand column. Fields
    default to 0 or empty string, except the strand which defaults to '+'.
    """

    __slots__ = ["_fields"]

    def __init__(self, line=None, _len=None):
        """Constructs a BED record from a line of text. The length of the
        object matches the number of columns in the input line; in the case
        incompatible values, a BEDError exception is raised.

        The len parameter is unused, and included only for compatibility with
        pysam parser objects, such as 'asBED'. No minimum number of columns are
        required, and it is possible to construct an empty bed record.
        """
        self._fields = []

        if line:
            line = line.rstrip("\r\n").split("\t")
            for column, (value, func) in enumerate(zip(line, _BED_TYPES)):
                try:
                    self._fields.append(func(value))
                except ValueError:
                    raise BEDError(
                        "Error parsing column %i in BED record "
                        "(%r); expected type %s, but found %r."
                        % (column, "\t".join(line), func.__name__, value)
                    )

            if len(line) > len(self._fields):
                self._fields.extend(line[len(self._fields) :])

    def freeze(self):
        record = BEDRecord()
        record._fields = tuple(self._fields)
        return record

    def __copy__(self):
        """Needed for copy.copy to work correctly as expected."""
        record = BEDRecord()
        record._fields = copy.copy(self._fields)
        return record

    def __len__(self):
        """Returns the number of fields in the record; 0 .. N."""
        return len(self._fields)

    def __str__(self):
        """Returns a string suitable for writing to a .bed file."""
        return "\t".join(str(value) for value in self._fields)

    def __repr__(self):
        """Returns a printable representation of the record."""
        fields = []
        for name, value in zip(_BED_KEYS, self._fields):
            fields.append("%s=%r" % (name, value))

        fields.extend(repr(value) for value in self._fields[len(_BED_KEYS) :])

        return "BEDRecord(%s)" % (", ".join(fields))

    def __getitem__(self, index):
        return self._fields[index]

    def __setitem__(self, index, value):
        if len(self._fields) <= index:
            defaults = _BED_DEFAULTS[len(self._fields) : index + 1]
            self._fields.extend(defaults)
            while len(self._fields) <= index:
                self._fields.append("")

        if index < len(_BED_TYPES):
            if type(_BED_TYPES[index]) is type:
                if not isinstance(value, _BED_TYPES[index]):
                    raise ValueError(
                        "Expected %s for BED field %i, got %r"
                        % (_BED_TYPES[index].__name__, index + 1, value)
                    )
            else:
                value = _BED_TYPES[index](value)

        self._fields[index] = value

    def __lt__(self, other):
        if not isinstance(other, BEDRecord):
            return NotImplemented

        return self._fields < other._fields

    @classmethod
    def _set_properties(cls):
        for index, name in enumerate(_BED_KEYS):
            setattr(cls, name, cls._new_attr(index))

    @classmethod
    def _new_attr(cls, index):
        """Returns an getter / setter property for the given value."""

        def _get(self):
            return self._fields[index]

        def _set(self, value):
            self[index] = value

        return property(_get, _set)


# Fill out properties for BEDRecord
BEDRecord._set_properties()


def read_bed_file(filename, min_columns=3, contigs=None):
    """Parses a (gzip/bzip2 compressed) BED file, and yields a sequence of
    records. Comments and empty lines are skipped. If the number of columns in
    the bed record is less than the specified ('min_columns'), a BEDError is
    raised. If a dictionary of {contig: length} is supplied, and min_columns
    is at least 6, then the coordinates are validated against the known contig
    lengths.
    """
    if min_columns < 3:
        raise ValueError("'min_columns' must be >= 3 in 'read_bed_file'")

    infinite = float("inf")
    handle = None
    try:
        handle = open_ro(filename)

        for (line_num, line) in enumerate(handle):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            try:
                bed = BEDRecord(line)
            except ValueError as error:
                raise BEDError(
                    "Error parsing line %i in regions file:\n"
                    "  Path = %r\n  Line = %r\n\n%s"
                    % (line_num + 1, filename, line, error)
                )

            if len(bed) < min_columns:
                url = "http://genome.ucsc.edu/FAQ/FAQformat.html#format1"
                name = repr(bed.name) if len(bed) > 3 else "unnamed record"
                raise BEDError(
                    "Region at line #%i (%s) does not "
                    "contain the expected number of fields; "
                    "the first %i fields are required. C.f. "
                    "defination at\n   %s\n\nPath = %r"
                    % (line_num, name, min_columns, url, filename)
                )

            if contigs is None:
                contig_len = infinite
            else:
                contig_len = contigs.get(bed.contig)

            if contig_len is None:
                raise BEDError(
                    "Regions file contains contig not found "
                    "in reference:\n  Path = %r\n  Contig = "
                    "%r\n\nPlease ensure that all contig "
                    "names match the reference names!" % (filename, bed.contig)
                )
            elif not (0 <= bed.start < bed.end <= contig_len):
                raise BEDError(
                    "Regions file contains invalid region:\n"
                    "  Path   = %r\n  Contig = %r\n"
                    "  Start  = %s\n  End    = %s\n\n"
                    "Expected 0 <= Start < End <= %i!"
                    % (filename, bed.contig, bed.start, bed.end, contig_len)
                )

            yield bed
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

    references = bamfile.references
    indices = dict(zip(references, range(len(references))))

    def _by_bam_layout(region):
        return (indices[region.contig], region.start, region.end)

    regions.sort(key=_by_bam_layout)


def pad_bed_records(records, padding, max_sizes={}):
    infinity = float("inf")
    for record in records:
        max_length = max_sizes.get(record.contig, infinity)
        record.start = max(0, record.start - padding)
        record.end = min(record.end + padding, max_length)


def merge_bed_records(records):
    records = sorted(records)
    if not records:
        return []

    last_record = BEDRecord()
    last_record._fields = records[0]._fields[:3]

    results = [last_record]
    for record in records:
        if last_record.contig != record.contig or last_record.end < record.start:
            last_record = BEDRecord()
            last_record._fields = record._fields[:3]
            results.append(last_record)
        else:
            last_record.end = record.end

    return results
