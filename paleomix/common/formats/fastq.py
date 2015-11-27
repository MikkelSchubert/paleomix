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


# Quality score offsets for Phred (or similar) scores in FASTQ reads (33 or 64)
OFFSET_33 = 33
OFFSET_64 = 64
# Quality score found in both ranges that are unique to each offset,
# suggesting that the list contains mixed quality offsets, or invalid data.
OFFSET_BOTH = "BOTH"
# No quality scores in the expected range of ASCII values (33 .. 105)
OFFSET_MISSING = "MISSING"
# Quality scores are in the ASCII range 59 .. 74, which could signify
# low-quality reads with offset 64, or high-quality reads with offset 33
OFFSET_AMBIGIOUS = "AMBIGIOUS"


def classify_quality_strings(quality_strings):
    """Takes a sequence of quality strings from FASTQ"""
    counts = [0] * 256
    for read in quality_strings:
        for char in read:
            counts[ord(char)] += 1

    return _identify_format(counts)


def _identify_format(counts):
    """Given a list representing counts of ASCII characters found in one or
    more FASTQ quality strins, this function attempts to identify the offset
    used to encode the quality scores.

    The following constants may be returned:
      - OFFSET_33: Offset identified as being 33
      - OFFSET_64: Offset identified as being 64
      - OFFSET_BOTH: Both offset 33 and 64 found, mixed file? (error)
      - OFFSET_MISSING: No quality scores found, wrong file? (error)
      - OFFSET_AMBIGIOUS: Qualities could be either offset. (warning)
    """
    # The range of scores that can unambigiously be identified
    # as belonging to Phred scores with offset 33 or 64. Scores
    # in between could potentially signify either offset
    # See e.g. http://en.wikipedia.org/wiki/FASTQ_format#Encoding
    has_offset_33_scores = any(counts[33:59])
    has_ambigious_scores = any(counts[59:75])
    has_offset_64_scores = any(counts[75:105])

    if has_offset_33_scores:
        if has_offset_64_scores:
            return OFFSET_BOTH
        return OFFSET_33
    elif has_offset_64_scores:
        return OFFSET_64
    elif has_ambigious_scores:
        return OFFSET_AMBIGIOUS
    return OFFSET_MISSING
