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
import sys
from itertools import groupby
from typing import TYPE_CHECKING, Iterable, Iterator

import pysam

from paleomix.common.argparse import ArgumentParser
from paleomix.common.text import padded_table
from paleomix.common.timer import BAMTimer

if TYPE_CHECKING:
    from pysam import AlignedSegment


def key_with_name(record: AlignedSegment) -> tuple[str | None, int, int]:
    return (record.query_name, record.reference_id, record.reference_start)


def key_without_name(record: AlignedSegment) -> tuple[int, int]:
    return (record.reference_id, record.reference_start)


def read_type(record: AlignedSegment) -> str:
    if record.is_paired:
        if record.is_read1:
            return "Mate 1 read in"
        elif record.is_read2:
            return "Mate 2 read in"
        else:
            return "Unpaired read in"
    else:
        return "Unpaired read in"


def print_duplicates(
    handle: pysam.AlignmentFile,
    duplicates: list[AlignedSegment],
) -> None:
    record = duplicates[0]
    print(
        f"Found {len(duplicates)} duplicates at "
        f"{record.reference_name}:{record.reference_start}"
    )

    samples = {}
    for rgroup in handle.header.get("RG", ()):
        samples[rgroup["ID"]] = (
            rgroup.get("SM", "?"),
            rgroup.get("LB", "?"),
            rgroup.get("PU", "?"),
        )

    groups: dict[tuple[str, str, str], list[AlignedSegment]] = collections.defaultdict(
        list
    )

    for record in duplicates:
        try:
            key = record.get_tag("RG")
        except KeyError:
            key = None

        key = samples.get(key, ("?", "?", "?"))
        groups[key].append(record)

    for (sample, library, run), records in sorted(groups.items()):
        rows: list[list[str | None]] = []

        for record in records:
            query_sequence = record.query_sequence
            if len(query_sequence) > 50:
                query_sequence = query_sequence[:50] + "..."

            query_qualities = "".join(chr(33 + qual) for qual in record.query_qualities)
            if len(query_qualities) > 50:
                query_qualities = query_qualities[:50] + "..."

            rows.append(
                [
                    read_type(record),
                    sample,
                    library,
                    run + ":",
                    record.query_name,
                    query_sequence,
                    query_qualities,
                ]
            )

        for line in padded_table(rows, min_padding=1):
            print(" ", line)

    print()


def complexity(record: AlignedSegment) -> int:
    last_nuc = None
    complexity = 0

    last_nuc = None
    for nuc in record.query_sequence:
        if nuc != last_nuc:
            complexity += 1
            last_nuc = nuc

    return complexity


def filter_records(records: Iterable[AlignedSegment]) -> Iterator[AlignedSegment]:
    for record in records:
        # Ignore unaligned / supplementary / secondary alignments
        if not record.flag & 0x904:
            yield record


def group_by_sequences(
    records: Iterable[AlignedSegment],
) -> Iterator[list[AlignedSegment]]:
    by_sequence: dict[tuple[bool, str | None], list[AlignedSegment]] = (
        collections.defaultdict(list)
    )
    for record in records:
        query_sequence = record.query_sequence
        if query_sequence:
            by_sequence[(record.is_reverse, record.query_sequence)].append(record)

    for group in by_sequence.values():
        if len(group) > 1:
            yield group


def group_by_qualities(
    records: Iterable[AlignedSegment],
) -> Iterator[list[AlignedSegment]]:
    by_qualities: dict[tuple[int, ...] | None, list[AlignedSegment]] = (
        collections.defaultdict(list)
    )
    for record in records:
        query_qualities = record.query_qualities
        if query_qualities is not None:
            query_qualities = tuple(query_qualities)

        by_qualities[query_qualities].append(record)

    for group in by_qualities.values():
        if len(group) > 1:
            yield group


def filter_candidate_duplicates(
    candidates: Iterable[AlignedSegment], min_complexity: float
) -> Iterator[tuple[list[AlignedSegment], list[AlignedSegment]]]:
    for group in group_by_sequences(filter_records(candidates)):
        for group in group_by_qualities(group):
            filtered_records = group
            if min_complexity > 0:
                filtered_records = [
                    record for record in group if complexity(record) >= min_complexity
                ]

            yield group, filtered_records


def build_parser() -> ArgumentParser:
    parser = ArgumentParser(
        prog="paleomix dupcheck",
        description="Attempt to detect reads included multiple times as input based "
        "on the presence of reads with identical names, sequences, and quality scores. "
        "Alternatively, tool can optionally can detect identical reads based only on "
        "the read sequence and quality scores. The input data is assumed to be mapped "
        "and coordinate sorted BAM alignments.",
    )

    parser.add_argument(
        "bam",
        nargs="?",
        help="BAM file to check for duplicated data [STDIN]",
    )
    parser.add_argument(
        "--max-output",
        type=int,
        default=100,
        help="Print as this many duplicated reads across all positions. Set to zero to "
        "print all duplicated reads",
    )
    parser.add_argument(
        "--min-complexity",
        type=float,
        default=0.0,
        help="Minimum complexity of reads required for them to be considered possible "
        "duplicate data. This is intended to exclude low-complexity and easily reads "
        "produced by empty spots and the like, which can easily result in identical "
        "sequences",
    )
    parser.add_argument(
        "--any-query-name",
        action="store_true",
        help="If set reads are can be considered duplicates despite their names "
        "differing. Sequences and quality scores must still be identical.",
    )

    return parser


def main(argv: list[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.bam is None:
        args.bam = "-"

        if sys.stdin.isatty():
            parser.print_help()
            return 0

    if not args.max_output:
        args.max_output = float("inf")

    key_function = key_without_name if args.any_query_name else key_with_name

    total_duplicates = 0
    total_filtered_duplicates = 0
    truncated_output_warning = args.max_output == 0
    with pysam.AlignmentFile(args.bam) as handle:
        for _, reads in groupby(BAMTimer(handle), key_function):
            filtered_duplicates = filter_candidate_duplicates(
                candidates=reads, min_complexity=args.min_complexity
            )

            for duplicates, high_complexity_duplicates in filtered_duplicates:
                total_duplicates += bool(duplicates)
                total_filtered_duplicates += bool(high_complexity_duplicates)
                if args.max_output > 0:
                    args.max_output -= 1

                    print_duplicates(handle, high_complexity_duplicates)
                elif not truncated_output_warning:
                    truncated_output_warning = True
                    print(
                        f"No further errors will be printed (max = {args.max_output}).."
                    )

    if not total_filtered_duplicates:
        print("No duplicate data found.")
        return 0

    print()
    print()
    print(
        f"A total of {total_duplicates} reads with duplicates were found, "
        f"with {total_filtered_duplicates} remaining "
    )
    print("after excluding low-complexity reads.")
    print()
    print("This indicates that the same data has been included multiple times in")
    print("the project. This can be because multiple copies of the same files were")
    print("used, or because one or more files contain multiple copies of the same ")
    print("reads.")

    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
