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
import sys

from paleomix.common.utilities import Immutable, TotallyOrdered
from paleomix.common.fileutils import open_ro
from paleomix.common.formats._common import FormatError


class FASTQError(FormatError):
    pass


class FASTQ(TotallyOrdered, Immutable):
    def __init__(self, name, meta, sequence, qualities):
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

    def write(self, fileobj=sys.stdout):
        """Writes a FASTQ record to fileobj."""
        name = self.name
        if self.meta:
            name = "%s %s" % (name, self.meta)

        fileobj.write("@%s\n%s\n+\n%s\n" % (name, self.sequence, self.qualities))

    @classmethod
    def from_lines(cls, lines):
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
                raise FASTQError("Invalid FASTQ header: %r" % (header,))

            header_fields = header.split(None, 1)
            name = header_fields[0]
            meta = header_fields[1] if len(header_fields) == 2 else None

            try:
                sequence = next(lines_iter).rstrip()
                separator = next(lines_iter).rstrip()
                qualities = next(lines_iter).rstrip()
            except StopIteration:
                raise FASTQError("Partial FASTQ record: %r" % (name,))

            if not separator.startswith("+"):
                raise FASTQError(
                    "Invalid FASTQ separator for %r; expected '+', found %r"
                    % (name, separator)
                )
            elif len(sequence) != len(qualities):
                raise FASTQError(
                    "Sequence length does not match qualities length for %r" % (name,)
                )

            yield FASTQ(
                name=name[1:], meta=meta, sequence=sequence, qualities=qualities
            )

    @classmethod
    def from_file(cls, filename):
        """Reads an unindexed FASTQ file, returning a sequence of
        tuples containing the name and sequence of each entry in
        the file. The FASTQ file may be GZIP/BZ2 compressed."""
        with open_ro(filename) as handle:
            yield from FASTQ.from_lines(handle)

    def __lt__(self, other):
        if not isinstance(other, FASTQ):
            return NotImplemented

        return (self.name, self.meta, self.sequence, self.qualities) < (
            other.name,
            other.meta,
            other.sequence,
            other.qualities,
        )

    def __hash__(self):
        return hash((self.name, self.meta, self.sequence, self.qualities))

    def __repr__(self):
        return "FASTQ(%r, %r, %r, %r)" % (
            self.name,
            self.meta,
            self.sequence,
            self.qualities,
        )


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

    def __init__(self):
        self._qualities = set()

    def update(self, record):
        self._qualities.update(record.qualities)

    def offsets(self):
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
                return FASTQualities.BOTH
            return FASTQualities.OFFSET_33
        elif has_offset_64_scores:
            return FASTQualities.OFFSET_64
        elif has_ambigious_scores:
            return FASTQualities.AMBIGIOUS

        return FASTQualities.MISSING
