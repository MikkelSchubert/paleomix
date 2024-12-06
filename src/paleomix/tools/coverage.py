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

import copy
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, Iterator, Mapping, Sequence, Union

from paleomix.common.bamfiles import BAMRegion, BAMRegionsIter
from paleomix.common.fileutils import file_or_stdout
from paleomix.common.timer import BAMTimer
from paleomix.tools.bam_stats.common import (
    Args,
    BAMStatsError,
    collect_readgroups,
    collect_references,
    main_wrapper,
)

if TYPE_CHECKING:
    import pysam
    from typing_extensions import TypeAlias

    CoverageTable: TypeAlias = defaultdict[
        str, defaultdict[str, defaultdict[str, Dict[str, "CoverageStats"]]]
    ]
    CoverageTableComponent: TypeAlias = Mapping[
        str, Union["CoverageStats", "CoverageTableComponent"]
    ]
    CountsTable: TypeAlias = Dict[str, Dict[str | None, "CoverageStats"]]


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


@dataclass
class CoverageStats:
    se: int = 0
    pe_1: int = 0
    pe_2: int = 0
    collapsed: int = 0
    hits: int = 0
    m: int = 0
    i: int = 0
    d: int = 0
    size: int = 0

    def add(self, other: CoverageStats) -> None:
        self.se += other.se
        self.pe_1 += other.pe_1
        self.pe_2 += other.pe_2
        self.collapsed += other.collapsed
        self.hits += other.hits
        self.m += other.m
        self.i += other.i
        self.d += other.d
        self.size += other.size

    def has_values(self) -> bool:
        return bool(
            self.se
            or self.pe_1
            or self.pe_2
            or self.collapsed
            or self.hits
            or self.m
            or self.i
            or self.d
        )


def _calculate_totals(table: CoverageTable) -> CoverageTable:
    lengths: dict[str, int] = {}
    for samples in table.values():
        for libraries in samples.values():
            for contigs in libraries.values():
                for name, contig in contigs.items():
                    size = lengths.get(name)
                    if (size is not None) and (size != contig.size):
                        raise BAMStatsError(name)
                    lengths[name] = contig.size

    for name, samples in sorted(table.items()):
        for sample, libraries in sorted(samples.items()):
            for library, contigs in sorted(libraries.items()):
                table[name][sample][library] = _calculate_totals_in(contigs, lengths)

            table[name][sample]["*"] = _calculate_totals_in(libraries, lengths)

        table[name]["*"]["*"] = _calculate_totals_in(table, lengths)
    return table


def _calculate_totals_in(
    tables: CoverageTableComponent,
    lengths: dict[str, int],
) -> dict[str, CoverageStats]:
    totals: dict[str, CoverageStats] = defaultdict(CoverageStats)
    total_size = sum(lengths.values())

    subtables = list(tables.items())
    while subtables:
        subtable_key, subtable = subtables.pop()
        if subtable_key == "*":
            totals[subtable_key].size = total_size
        elif isinstance(subtable, CoverageStats):
            totals[subtable_key].add(subtable)
            totals["*"].add(subtable)

            totals[subtable_key].size = lengths[subtable_key]
            totals["*"].size = total_size
        else:
            subtables.extend(subtable.items())

    return dict(totals)


def build_rows(table: CoverageTable) -> Iterator[Sequence[str | int | float]]:
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
                        subtable.size,
                        subtable.se
                        + subtable.pe_1
                        + subtable.pe_2
                        + subtable.collapsed,
                        subtable.se,
                        subtable.pe_1,
                        subtable.pe_2,
                        subtable.collapsed,
                        subtable.m,
                        subtable.i,
                        subtable.d,
                        float(subtable.m) / subtable.size,
                    ]
                    yield row
                yield "#"
            yield "#"


def write_table(table: CoverageTable, filename: str) -> None:
    table = _calculate_totals(table)
    rows = build_rows(table)

    with file_or_stdout(filename) as out:
        print(TABLE_HEADER, end="", file=out)
        for line in rows:
            print(*line, sep="\t", file=out)


def build_region_template(
    args: Args,
    handle: pysam.AlignmentFile,
) -> dict[str | None, CoverageStats]:
    template: dict[str | None, CoverageStats] = {}
    for readgroup in collect_readgroups(args, handle):
        template[readgroup.key] = CoverageStats()
    return template


def get_region_table(
    counts: CountsTable,
    region: str,
    template: dict[str | None, CoverageStats],
) -> dict[str | None, CoverageStats]:
    subtable = counts.get(region)
    if subtable is not None:
        return subtable

    subtable = copy.deepcopy(template)
    counts[region] = subtable

    return subtable


def build_table(
    args: Args,
    handle: pysam.AlignmentFile,
    counts: CountsTable,
) -> CoverageTable:
    references = collect_references(args, handle)

    table: CoverageTable = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for readgroup in collect_readgroups(args, handle):
        # Exclude counts for reads with no read-groups, if none such were seen
        if readgroup.key is None and not args.ignore_readgroups:
            for reference in references:
                stats = counts.get(reference, {}).get(readgroup.key)
                if stats is not None and stats.has_values():
                    break
            else:
                continue

        for reference, size in references.items():
            subtable = table[args.target_name][readgroup.sample][readgroup.library]
            try:
                stats = subtable[reference]
            except KeyError:
                stats = subtable[reference] = CoverageStats(size=size)

            statistics = counts.get(reference, {}).get(readgroup.key)
            if statistics is not None:
                stats.add(statistics)

    return table


def print_table(
    args: Args,
    handle: pysam.AlignmentFile,
    counts: CountsTable,
) -> None:
    table = build_table(args, handle, counts)
    write_table(table, args.outfile)


##############################################################################
##############################################################################


def process_record(
    subtable: CoverageStats,
    record: pysam.AlignedSegment,
    flags: int,
    region: BAMRegion,
) -> None:
    qname = record.query_name
    if qname and qname.startswith(("M_", "MT_")):
        subtable.collapsed += 1
    elif flags & 0x40:  # first of pair
        subtable.pe_1 += 1
    elif flags & 0x80:  # second of pair
        subtable.pe_2 += 1
    else:  # Singleton
        subtable.se += 1

    position = record.query_alignment_start
    start = region.start

    end = region.end
    cigartuples = record.cigartuples
    # Reads are filterd and should normally have cigars/end-coordinates
    if cigartuples is not None and end is not None:
        for cigar, num in cigartuples:
            left = min(max(position, start), end)
            right = min(max(position + num, start), end)
            bases_in_region = right - left
            assert 0 <= bases_in_region <= num

            # 0 = 'M', 1 = 'I', 2 = 'D', 7 = '=', 8 = 'X'
            if cigar in (0, 7, 8):  # M, =, X
                subtable.m += bases_in_region
                position += num
            elif cigar == 1:  # I
                subtable.i += bases_in_region
            elif cigar == 2:  # D
                subtable.d += bases_in_region
                position += num
            elif cigar == 3:  # N
                position += num


def process_file(args: Args, handle: pysam.AlignmentFile) -> int:
    timer = BAMTimer(handle)

    counts: CountsTable = {}
    last_tid = 0
    region_template = build_region_template(args, handle)
    for region in BAMRegionsIter(handle, args.regions):
        if region.name is None:
            # Trailing unmapped reads
            break

        name = region.name
        if not args.regions and (handle.nreferences > args.max_contigs):
            name = "<Genome>"

        last_pos = 0
        region_table = get_region_table(counts, name, region_template)
        for position, records in region:
            for record in records:
                readgroup = args.get_readgroup_func(record)
                readgroup_table = region_table.get(readgroup)
                if readgroup_table is None:
                    # Unknown readgroups are treated as missing readgroups
                    readgroup_table = region_table[None]

                process_record(readgroup_table, record, record.flag, region)
                timer.increment()

            if (region.tid, position) < (last_tid, last_pos):
                sys.stderr.write("ERROR: Input BAM file is unsorted\n")
                return 1

            last_pos = position
            last_tid = region.tid

    timer.finalize()

    print_table(args, handle, counts)

    return 0


def main(argv: list[str]) -> int:
    return main_wrapper(process_file, argv, ".coverage")


##############################################################################
##############################################################################

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
