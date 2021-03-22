#!/usr/bin/env python3
import json
import logging
import sys

from collections import defaultdict
from pathlib import Path

import pysam

from paleomix.common.argparse import ArgumentParser
from paleomix.common.logging import initialize_console_logging
from paleomix.common.bamfiles import (
    BAM_SEGMENTED,
    BAM_SUPPLEMENTARY_ALIGNMENT,
    BAM_PCR_DUPLICATE,
    BAM_QUALITY_CONTROL_FAILED,
    BAM_SECONDARY_ALIGNMENT,
    BAM_NEXT_IS_UNMAPPED,
    BAM_READ_IS_UNMAPPED,
    BAM_PROPER_SEGMENTS,
)
from paleomix.common.timer import BAMTimer

_SE_FILTER_AND_NAMES = (
    ("unmapped", BAM_READ_IS_UNMAPPED),
    ("secondary", BAM_SECONDARY_ALIGNMENT),
    ("qc_failed", BAM_QUALITY_CONTROL_FAILED),
    ("pcr_duplicate", BAM_PCR_DUPLICATE),
    ("supplementary", BAM_SUPPLEMENTARY_ALIGNMENT),
)

_PE_FILTER_AND_NAMES = _SE_FILTER_AND_NAMES + (
    ("unmapped_mate", BAM_NEXT_IS_UNMAPPED),
    # Proper segments is a good thing, so we flip the bit below to catch improper pairs
    ("improper_pair", BAM_PROPER_SEGMENTS),
)

_SE_FILTERS = sum(value for _, value in _SE_FILTER_AND_NAMES)
_PE_FILTERS = sum(value for _, value in _PE_FILTER_AND_NAMES)


def is_proper_alignment(flag):
    if flag & BAM_SEGMENTED:
        return not flag & _PE_FILTERS
    else:
        return not flag & _SE_FILTERS


def calculate_flag_statistics(statistics, readgroups, flag_counts):
    for readgroup, counts in flag_counts.items():
        description = readgroups[readgroup]
        stats = statistics[description]

        for flag, count in counts.items():
            if flag & BAM_SEGMENTED:
                masks = _PE_FILTER_AND_NAMES
            else:
                masks = _SE_FILTER_AND_NAMES

            passes_all = True
            stats["total"] += count
            for label, mask in masks:
                if flag & mask:
                    stats[label] += count
                    passes_all = False

            if passes_all:
                stats["total_passed"] += count
            else:
                stats["total_failed"] += count


def calculate_cigar_statistics(statistics, readgroups, cigar_counts, genome_size, key):
    for readgroup, counts in cigar_counts.items():
        total_matches = 0
        for cigar, count in counts.items():
            for cigar_op, num in cigar:
                if cigar_op in (0, 7, 8):
                    total_matches += num

        description = readgroups[readgroup]
        stats = statistics[description]
        stats[key + "bases"] = total_matches
        stats[key + "coverage"] = total_matches / genome_size


def calculate_statistics(handle, flag_counts, passed_cigar_counts, failed_cigar_counts):
    readgroups = {}
    for readgroup in handle.header["RG"]:
        readgroups[readgroup["ID"]] = readgroup["DS"]

    genome_size = sum(handle.lengths)

    statistics = defaultdict(lambda: defaultdict(int))
    calculate_flag_statistics(statistics, readgroups, flag_counts)
    calculate_cigar_statistics(
        statistics=statistics,
        readgroups=readgroups,
        cigar_counts=passed_cigar_counts,
        genome_size=genome_size,
        key="total_passed_",
    )
    calculate_cigar_statistics(
        statistics=statistics,
        readgroups=readgroups,
        cigar_counts=failed_cigar_counts,
        genome_size=genome_size,
        key="total_failed_",
    )

    totals = defaultdict(int)
    for counts in tuple(statistics.values()):
        for key, value in counts.items():
            totals[key] += value

    statistics["*"] = totals

    return statistics


def parse_args(argv):
    parser = ArgumentParser("paleomix ngs:finalize_bam")
    parser.add_argument("in_bam", type=Path)
    parser.add_argument("--out-passed", type=Path)
    parser.add_argument("--out-failed", type=Path)
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--threads", default=1, type=int)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    initialize_console_logging()

    log = logging.getLogger("ngs:finalize_bam")
    if not (args.out_passed or args.out_failed or args.out_json):
        log.error("No --out-* arguments; please specify at least one output file")
        return 1

    log.info("Reading BAM records from %r", args.in_bam)

    # Counts of flags per read group
    if args.out_json:
        flag_counts = defaultdict(lambda: defaultdict(int))
        passed_cigar_counts = defaultdict(lambda: defaultdict(int))
        failed_cigar_counts = defaultdict(lambda: defaultdict(int))

    log.info("Reading alignents from %s", args.in_bam)
    in_bam = pysam.AlignmentFile(args.in_bam, threads=args.threads)

    if args.out_passed:
        log.info("Writing proper alignents to %r", args.out_failed)
        out_passed = pysam.AlignmentFile(
            args.out_passed, "wb", template=in_bam, threads=args.threads
        )
    else:
        out_passed = None

    if args.out_failed:
        log.info("Writing failed alignents/unaligned reads to %r", args.out_failed)
        out_failed = pysam.AlignmentFile(
            args.out_failed, "wb", template=in_bam, threads=args.threads
        )
    else:
        out_failed = None

    for record in BAMTimer(in_bam):
        # Flip bit so that 1 represents improper segments
        flag = record.flag ^ BAM_PROPER_SEGMENTS
        if args.out_json:
            read_group = record.get_tag("RG")
            flag_counts[read_group][flag] += 1

        if is_proper_alignment(flag):
            if out_passed:
                out_passed.write(record)

            if args.out_json:
                passed_cigar_counts[read_group][tuple(record.cigartuples)] += 1
        else:
            if out_failed:
                out_failed.write(record)

            if args.out_json:
                cigar = record.cigartuples
                if cigar is not None:
                    failed_cigar_counts[read_group][tuple(record.cigartuples)] += 1

    if out_passed:
        out_passed.close()

    if out_failed:
        out_failed.closer()

    if args.out_json:
        lengths = in_bam.lengths
        statistics = calculate_statistics(
            handle=in_bam,
            flag_counts=flag_counts,
            passed_cigar_counts=passed_cigar_counts,
            failed_cigar_counts=failed_cigar_counts,
        )

        with args.out_json.open("wt") as handle:
            json.dump(
                {
                    "input": str(args.in_bam),
                    "output_passed": str(args.out_passed),
                    "output_failed": str(args.out_failed),
                    "statistics": statistics,
                    "genome": {
                        "ncontigs": len(lengths),
                        "size": sum(lengths),
                    },
                },
                handle,
            )

    in_bam.close()

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv[1:]))
