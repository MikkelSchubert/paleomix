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
import collections
import sys

from paleomix.common.text import parse_padded_table
from paleomix.common.utilities import get_in, set_in
from paleomix.tools.bam_stats.common import BAMStatsError

##############################################################################
##############################################################################
##


class ReadGroup:
    __slots__ = ["SE", "PE_1", "PE_2", "Collapsed", "Hits", "M", "I", "D", "Size"]

    def __init__(self):
        self.SE = 0
        self.PE_1 = 0
        self.PE_2 = 0
        self.Collapsed = 0
        self.Hits = 0
        self.M = 0
        self.I = 0
        self.D = 0
        self.Size = 0

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def add(self, other):
        self.SE += other.SE
        self.PE_1 += other.PE_1
        self.PE_2 += other.PE_2
        self.Collapsed += other.Collapsed
        self.Hits += other.Hits
        self.M += other.M
        self.I += other.I
        self.D += other.D
        self.Size += other.Size

    def has_values(self):
        for slot in ReadGroup.__slots__:
            if slot != "Size" and getattr(self, slot):
                return True

        return False


# Header prepended to output tables
TABLE_HEADER = """# Columns:
#   Contig:    Contig, chromosome, or feature for which a depth histogram was
#              created. Unnamed features are named after the chromosome or
#              contig on which they are located, with a star appended. For
#              example "chr1*".
#   Size:      The total size of the region. Multiple features with the same
#              name are combined into one row, with the size representing to
#              total of these. Note that overlapping bases are counted 2 (or
#              more) times.
#   Hits:      Sum of SE, PE_1, and PE_2 hits. Note that supplementary
#              alignments, duplicates, reads that failed QC, secondary
#              alignments, and unmapped reads are ignored.
#   SE, PE_*:  Number of Single Ended, and Pair Ended (mate 1 and 2) hits
#              overlapping the current contig or intervals. Note that a hit
#              may be counted multiple times if it overlaps multiple intervals
#   Collapsed: Number of hits for PE pair collapsed into a single read.
#   M, I, D:   Number of aligned (M), inserted (I) and deleted (D) bases
#              relative to references.
#   Coverage:  Average number of bases covering each position in the
#              contig(s)/intervals(s).
"""


def calculate_totals(table):
    lengths = {}
    for samples in table.values():
        for libraries in samples.values():
            for contigs in libraries.values():
                for name, contig in contigs.items():
                    size = lengths.get(name)
                    if (size is not None) and (size != contig.Size):
                        raise BAMStatsError(name)
                    lengths[name] = contig.Size

    for name, samples in sorted(table.items()):
        for sample, libraries in sorted(samples.items()):
            for library, contigs in sorted(libraries.items()):
                totals = _calculate_totals_in(contigs, lengths)
                set_in(table, (name, sample, library), totals)

            totals = _calculate_totals_in(libraries, lengths)
            set_in(table, (name, sample, "*"), totals)

        set_in(table, (name, "*", "*"), _calculate_totals_in(table, lengths))
    return table


def build_rows(table):
    yield (
        "Name",
        "Sample",
        "Library",
        "Contig",
        "Size",
        "Hits",
        "SE",
        "PE_1",
        "PE_2",
        "Collapsed",
        "M",
        "I",
        "D",
        "Coverage",
    )

    for name, samples in sorted(table.items()):
        for sample, libraries in sorted(samples.items()):
            for library, contigs in sorted(libraries.items()):
                for contig, subtable in sorted(contigs.items()):
                    row = [
                        name,
                        sample,
                        library,
                        contig,
                        subtable.Size,
                        subtable.SE
                        + subtable.PE_1
                        + subtable.PE_2
                        + subtable.Collapsed,
                        subtable.SE,
                        subtable.PE_1,
                        subtable.PE_2,
                        subtable.Collapsed,
                        subtable.M,
                        subtable.I,
                        subtable.D,
                        float(subtable.M) / subtable.Size,
                    ]
                    yield row
                yield "#"
            yield "#"


def read_table(table, filename):
    with open(filename) as table_file:
        for record in parse_padded_table(table_file):
            key = (
                record["Name"],
                record["Sample"],
                record["Library"],
                record["Contig"],
            )
            if "*" in key:
                continue

            subtable = get_in(table, key)
            if subtable is None:
                subtable = ReadGroup()
                subtable.Size = int(record["Size"])
                set_in(table, key, subtable)

            assert int(subtable.Size) == int(record["Size"])
            for key in ReadGroup.__slots__:
                if key != "Size":
                    subtable[key] += int(record.get(key, 0))


def write_table(table, filename):
    table = calculate_totals(table)
    rows = build_rows(table)

    if filename == "-":
        output_handle = sys.stdout
    else:
        output_handle = open(filename, "w")

    try:
        output_handle.write(TABLE_HEADER)
        for line in rows:
            output_handle.write("\t".join(map(str, line)))
            output_handle.write("\n")
    finally:
        if output_handle is not sys.stdout:
            output_handle.close()


def _calculate_totals_in(tables, lengths):
    totals = collections.defaultdict(ReadGroup)
    total_size = sum(lengths.values())

    subtables = list(tables.items())
    while subtables:
        subtable_key, subtable = subtables.pop()
        if subtable_key == "*":
            totals[subtable_key].Size = total_size
        elif isinstance(subtable, ReadGroup):
            totals[subtable_key].add(subtable)
            totals["*"].add(subtable)

            totals[subtable_key].Size = lengths[subtable_key]
            totals["*"].Size = total_size
        else:
            subtables.extend(subtable.items())

    return dict(totals)
