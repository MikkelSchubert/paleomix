#!/usr/bin/env python3
# -*- coding: utf8 -*-
import collections
import sys

from itertools import groupby

import pysam

from paleomix.common.argparse import ArgumentParser
from paleomix.common.timer import BAMTimer
from paleomix.common.text import padded_table


def key_with_name(record):
    return (record.query_name, record.reference_id, record.reference_start)


def key_without_name(record):
    return (record.reference_id, record.reference_start)


def read_type(record):
    if record.is_paired:
        if record.is_read1:
            return "Mate 1 read in"
        elif record.is_read2:
            return "Mate 2 read in"
        else:
            return "Unpaired read in"
    else:
        return "Unpaired read in"


def print_duplicates(handle, duplicates):
    record = duplicates[0]
    print(
        "Found {} duplicates at {}:{}".format(
            len(duplicates),
            record.reference_name,
            record.reference_start,
        )
    )

    samples = {}
    for rgroup in handle.header.get("RG", ()):
        samples[rgroup["ID"]] = (
            rgroup.get("SM", "?"),
            rgroup.get("LB", "?"),
            rgroup.get("PU", "?"),
        )

    groups = collections.defaultdict(list)
    for record in duplicates:
        try:
            key = record.get_tag("RG")
        except KeyError:
            key = None

        key = samples.get(key, ("?", "?", "?"))
        groups[key].append(record)

    for (sample, library, run), records in sorted(groups.items()):
        rows = []

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


def complexity(record):
    last_nuc = None
    complexity = 0

    last_nuc = None
    for nuc in record.query_sequence:
        if nuc != last_nuc:
            complexity += 1
            last_nuc = nuc

    return complexity


def filter_records(records):
    for record in records:
        # Ignore unaligned / supplementary / secondary alignments
        if not record.flag & 0x904:
            yield record


def group_by_sequences(records):
    by_sequence = collections.defaultdict(list)
    for record in records:
        query_sequence = record.query_sequence
        if query_sequence:
            by_sequence[(record.is_reverse, record.query_sequence)].append(record)

    for _, group in by_sequence.items():
        if len(group) > 1:
            yield group


def group_by_qualities(records):
    by_qualities = collections.defaultdict(list)
    for record in records:
        query_qualities = record.query_qualities
        if query_qualities is not None:
            query_qualities = tuple(record.query_qualities)

        by_qualities[tuple(record.query_qualities)].append(record)

    for group in by_qualities.values():
        if len(group) > 1:
            yield group


def filter_candidate_duplicates(candidates, min_complexity):
    for group in group_by_sequences(filter_records(candidates)):
        for group in group_by_qualities(group):
            filtered_records = group
            if min_complexity > 0:
                filtered_records = [
                    record for record in group if complexity(record) >= min_complexity
                ]

            yield group, filtered_records


def parse_args():
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


def main(argv):
    parser = parse_args()
    args = parser.parse_args(argv)
    if args.bam is None:
        args.bam = "-"

        if sys.stdin.isatty():
            parser.print_help()
            return 0

    if not args.max_output:
        args.max_output = float("inf")

    if args.any_query_name:
        key_function = key_without_name
    else:
        key_function = key_with_name

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
                        "No further errors will be printed (max = {}) ..".format(
                            args.max_output
                        )
                    )

    if not total_filtered_duplicates:
        print("No duplicate data found.")
        return 0

    print()
    print()
    print(
        "A total of {} reads with duplicates were found, with {} remaining ".format(
            total_duplicates,
            total_filtered_duplicates,
        )
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
