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

import collections
import itertools
import sys
from collections.abc import Iterable, Iterator
from typing import NewType, Optional

import pysam
from typing_extensions import TypeAlias

from paleomix.common.bamfiles import BAMRegion, BAMRegionsIter
from paleomix.common.fileutils import file_or_stdout
from paleomix.common.timer import BAMTimer
from paleomix.tools.bam_stats.common import (
    Args,
    collect_readgroups,
    collect_references,
    main_wrapper,
)

SampleLibraryKey = NewType("SampleLibraryKey", tuple[str, str])
SampleLibraryID = NewType("SampleLibraryID", int)

CountsCache: TypeAlias = collections.deque[list[int]]
CountsDict: TypeAlias = dict[int, int]
TotalsKey: TypeAlias = tuple[str, str, str]
TotalsDict: TypeAlias = dict[TotalsKey, CountsDict]

ReadGroupToID: TypeAlias = dict[Optional[str], SampleLibraryID]
IDToLibrary: TypeAlias = list[SampleLibraryKey]


##############################################################################
##############################################################################
##

# Maximum depth to record, and hence the number of columns in output table
_MAX_DEPTH = 200
# Maximum number of count patterns (numbers of bases per library for a given
# site) to cache for bulk processing; see MappingsToTotals for implementation
_MAX_CACHE_SIZE = 10000


# Header prepended to output tables
_HEADER = """# Columns:
#   Contig:   Contig, chromosome, or feature for which a depth histogram was
#             created. Unnamed features are named after the chromosome or
#             contig on which they are located, with a star appended. For
#             example "chr1*". If the maximum number of contigs was exceeded,
#             these are collapsed into one meta-contig named "<Genome>".
#   Size:     The total size of the region. Multiple features with the same
#             name are combined into one row, with the size representing the
#             total for these. Note that overlapping bases are counted 2 (or
#             more) times.
#   MaxDepth: Maximum depth to use when calling SNPs, in order to exclude
#             (at least) the 0.5%% most extreme sites based on read depth,
#             not including sites with depth 0.
#   MD_*:     Fraction of sites with a minimum depth of 1-200.
#"""


##############################################################################
##############################################################################


class MappingToTotals:
    def __init__(
        self,
        totals: TotalsDict,
        region: BAMRegion,
        id_to_library: IDToLibrary,
    ) -> None:
        assert region.name is not None, "regions should have been named"

        self._region = region
        self._map_by_smlbid, self._totals_src_and_dst = self._build_mappings(
            totals, region.name, id_to_library
        )
        self._cache: dict[tuple[int, ...], int] = collections.defaultdict(int)

    def process_counts(
        self,
        counts: CountsCache,
        last_pos: int,
        cur_pos: float,
    ) -> None:
        start = self._region.start
        end = self._region.end
        if end is None:
            end = float("inf")

        # Pileups tends to contain identical stretches, so
        # try to avoid repeated lookups by aggregating these
        repeats = 1
        last_count: list[int] | None = None
        while counts and (last_pos < cur_pos):
            count = counts.popleft()
            if start <= last_pos < end:
                if count == last_count:
                    repeats += 1
                else:
                    if last_count is not None:
                        self._cache[tuple(last_count)] += repeats
                    last_count = count
                    repeats = 1
            last_pos += 1

        if last_count is not None:
            self._cache[tuple(last_count)] += repeats

        if len(self._cache) > _MAX_CACHE_SIZE:
            self.finalize()

    def finalize(self) -> None:
        """Process cached counts."""
        for counts, multiplier in self._cache.items():
            for smlbid, count in enumerate(counts):
                if count:
                    for lst in self._map_by_smlbid[smlbid]:
                        lst[0] += count

            for dst_counts, src_count in self._totals_src_and_dst:
                if src_count[0]:
                    dst_counts[src_count[0]] += multiplier
                    src_count[0] = 0

        self._cache.clear()

    @classmethod
    def _build_mappings(
        cls,
        totals: TotalsDict,
        name: str,
        id_to_library: IDToLibrary,
    ) -> tuple[list[tuple[CountsDict, ...]], list[tuple[CountsDict, CountsDict]]]:
        # Accumulators mapped by sample+library IDs
        totals_by_smlbid: list[tuple[CountsDict, ...]] = []
        # Accumulators mapped by the corresponding table keys
        totals_by_table_key: TotalsDict = {}

        for sm_key, lb_key in id_to_library:
            keys = [
                ("*", "*", "*"),
                (sm_key, "*", "*"),
                (sm_key, "*", name),
                (sm_key, lb_key, "*"),
                (sm_key, lb_key, name),
            ]

            mappings = cls._nonoverlapping_mappings(keys, totals, totals_by_table_key)
            totals_by_smlbid.append(mappings)

        totals_src_and_dst: list[tuple[CountsDict, CountsDict]] = []
        for key, dst in totals_by_table_key.items():
            totals_src_and_dst.append((totals[key], dst))

        return totals_by_smlbid, totals_src_and_dst

    @staticmethod
    def _nonoverlapping_mappings(
        keys: Iterable[TotalsKey],
        totals: TotalsDict,
        totals_by_table_key: TotalsDict,
    ) -> tuple[CountsDict, ...]:
        """Returns a tuple of accumulators for a given set of table keys. As
        multiple table keys may share the same accumulator (e.g. if there is
        only one sample, then sample "*" and that sample will be identical),
        the tuple of accumulators may contain fewer items than keys."""

        mapping: list[CountsDict] = []
        totals_used: set[int] = set()
        for key in keys:
            # Check that accumulator is not already included
            totals_id = id(totals[key])
            if totals_id not in totals_used:
                totals_used.add(totals_id)
                accumulator = totals_by_table_key.setdefault(key, {0: 0})
                mapping.append(accumulator)
        return tuple(mapping)


##############################################################################
##############################################################################


def calc_max_depth(counts: CountsDict) -> int | str:
    counts = dict(counts)
    counts.pop(0, None)

    running_total = sum(counts.values())
    if not running_total:
        return "NA"

    total = float(running_total)
    for index, count in sorted(counts.items()):
        # Stop when less than the 0.5% most extreme values are included
        if running_total / total < 0.005:
            # The max is inclusive, so return the depth just before this one
            return index - 1
        running_total -= count

    return "NA"


def print_table(args: Args, handle: pysam.AlignmentFile, totals: TotalsDict) -> None:
    lengths = collect_references(args, handle)

    with file_or_stdout(args.outfile) as out:
        print(_HEADER, file=out)
        for line in build_table(args.target_name, totals, lengths):
            print(*line, sep="\t", file=out)


def calculate_depth_pc(counts: CountsDict, length: float) -> Iterator[str]:
    final_counts = [0] * (_MAX_DEPTH + 1)
    for depth, count in counts.items():
        final_counts[min(_MAX_DEPTH, depth)] += count

    running_total = sum(final_counts)
    for count in final_counts[1:]:
        yield f"{running_total / length:.4f}"
        running_total -= count


def build_table(
    name: str,
    totals: TotalsDict,
    lengths: dict[str, int],
) -> Iterator[list[str] | str]:
    header = ["Name", "Sample", "Library", "Contig", "Size", "MaxDepth"]
    for index in range(1, _MAX_DEPTH + 1):
        header.append(f"MD_{index:03d}")

    yield header
    last_sm = last_lb = None
    for (sm_key, lb_key, ct_key), counts in sorted(totals.items()):
        if (sm_key != last_sm) and (last_sm is not None):
            yield "#"
            yield "#"
        elif (lb_key != last_lb) and (last_lb is not None):
            yield "#"
        last_sm, last_lb = sm_key, lb_key

        length = sum(lengths.values()) if ct_key == "*" else lengths[ct_key]

        row = [name, sm_key, lb_key, ct_key, str(length), str(calc_max_depth(counts))]
        row.extend(calculate_depth_pc(counts, length))

        yield row


##############################################################################
##############################################################################


def build_key_struct(
    args: Args, handle: pysam.AlignmentFile
) -> collections.defaultdict[str, set[str]]:
    structure: dict[str, set[str]] = collections.defaultdict(set)
    for readgroup in collect_readgroups(args, handle):
        structure[readgroup.sample].add(readgroup.library)

    return structure


def build_new_dicts(
    totals: TotalsDict,
    sample: str,
    library: str,
    contigs: Iterable[str],
) -> None:
    totals[(sample, library, "*")] = collections.defaultdict(int)
    for contig in contigs:
        totals[(sample, library, contig)] = collections.defaultdict(int)


def reuse_dicts(
    totals: TotalsDict,
    dst_sm: str,
    dst_lb: str,
    src_sm: str,
    src_lb: str,
    contigs: Iterable[str],
) -> None:
    totals[(dst_sm, dst_lb, "*")] = totals[(src_sm, src_lb, "*")]
    for contig in contigs:
        totals[(dst_sm, dst_lb, contig)] = totals[(src_sm, src_lb, contig)]


def build_totals_dict(args: Args, handle: pysam.AlignmentFile) -> TotalsDict:
    references = tuple(collect_references(args, handle))
    structure = build_key_struct(args, handle)

    totals: TotalsDict = {}
    for sm_key, libraries in structure.items():
        for lb_key in libraries:
            if len(references) == 1:
                key = references[0]
                counts: CountsDict = collections.defaultdict(int)
                totals[(sm_key, lb_key, key)] = counts
                totals[(sm_key, lb_key, "*")] = counts
            else:
                build_new_dicts(totals, sm_key, lb_key, references)

        if len(libraries) == 1:
            (key,) = libraries
            reuse_dicts(totals, sm_key, "*", sm_key, key, references)
        else:
            build_new_dicts(totals, sm_key, "*", references)

    if len(structure) == 1:
        (key,) = structure
        reuse_dicts(totals, "*", "*", key, "*", references)
    else:
        build_new_dicts(totals, "*", "*", references)

    return totals


def count_bases(
    args: Args,
    counts: CountsCache,
    record: pysam.AlignedSegment,
    rg_to_smlbid: ReadGroupToID,
    template: list[int],
) -> None:
    cigartuples = record.cigartuples
    if cigartuples is not None and record.reference_length is not None:
        for _ in range(record.reference_length - len(counts)):
            counts.append(list(template))

        key = rg_to_smlbid.get(args.get_readgroup_func(record))
        if key is None:
            # Unknown readgroups are treated as missing readgroups
            key = rg_to_smlbid[None]

        index = 0
        for cigar, count in cigartuples:
            if cigar in (0, 7, 8):
                for counter in itertools.islice(counts, index, index + count):
                    counter[key] += 1
                index += count
            elif cigar in (2, 3, 6):
                index += count


def build_rg_to_smlbid_keys(
    args: Args,
    handle: pysam.AlignmentFile,
) -> tuple[ReadGroupToID, IDToLibrary]:
    """Returns a dictionary which maps a readgroup ID to an index value,
    as well as a list containing a tuple (sample, library) corresponding
    to each index. Typically, this list will be shorter than the map of read-
    groups, as multiple read-groups will map to the same sample / library.
    """
    readgroup_to_id: dict[str | None, SampleLibraryID] = {}
    library_to_id: dict[SampleLibraryKey, SampleLibraryID] = {}
    id_to_library: list[SampleLibraryKey] = []
    for readgroup in collect_readgroups(args, handle):
        key = SampleLibraryKey((readgroup.sample, readgroup.library))
        if key not in library_to_id:
            library_to_id[key] = SampleLibraryID(len(library_to_id))
            id_to_library.append(key)

        readgroup_to_id[readgroup.key] = library_to_id[key]

    return readgroup_to_id, id_to_library


def process_file(args: Args, handle: pysam.AlignmentFile) -> int:
    timer = BAMTimer(handle)

    last_tid = 0
    totals = build_totals_dict(args, handle)
    rg_to_smlbid, id_to_library = build_rg_to_smlbid_keys(args, handle)
    template = [0] * len(id_to_library)

    for region in BAMRegionsIter(handle, args.regions):
        if region.name is None:
            # Trailing unmapped reads
            break
        elif not args.regions and (handle.nreferences > args.max_contigs):
            region.name = "<Genome>"

        last_pos = 0
        counts: CountsCache = collections.deque()
        mapping = MappingToTotals(totals, region, id_to_library)
        for position, records in region:
            mapping.process_counts(counts, last_pos, position)

            for record in records:
                timer.increment()
                count_bases(args, counts, record, rg_to_smlbid, template)

            if (region.tid, position) < (last_tid, last_pos):
                sys.stderr.write("ERROR: Input BAM file is unsorted\n")
                return 1

            last_pos = position
            last_tid = region.tid

        # Process columns in region after last read
        mapping.process_counts(counts, last_pos, float("inf"))
        mapping.finalize()
    timer.finalize()

    if not args.ignore_readgroups:
        # Exclude counts for reads with no read-groups, if none such were seen
        for (key, _, _), value in totals.items():
            if key == "<NA>" and value:
                break
        else:
            for key in list(totals):
                if key[0] == "<NA>":
                    totals.pop(key)

    print_table(args, handle, totals)

    return 0


def main(argv: list[str]) -> int:
    return main_wrapper(process_file, argv, ".depths")


##############################################################################
##############################################################################

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
