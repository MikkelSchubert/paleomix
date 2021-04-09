#!/usr/bin/env python3
import functools
import json
import logging
import sys

from collections import defaultdict
from itertools import groupby
from pathlib import Path

import pysam

from paleomix.common.argparse import ArgumentParser
from paleomix.common.logging import initialize_console_logging
from paleomix.common.bamfiles import (
    BAM_IS_FIRST_SEGMENT,
    BAM_IS_LAST_SEGMENT,
    BAM_NEXT_IS_UNMAPPED,
    BAM_PCR_DUPLICATE,
    BAM_PROPER_SEGMENTS,
    BAM_QUALITY_CONTROL_FAILED,
    BAM_READ_IS_UNMAPPED,
    BAM_SECONDARY_ALIGNMENT,
    BAM_SEGMENTED,
    BAM_SUPPLEMENTARY_ALIGNMENT,
)
from paleomix.common.timer import BAMTimer


def identify_read(record):
    return record.get_tag("RG")


@functools.lru_cache()
def count_mapped_bases(cigar):
    total_matches = 0
    for cigar_op, num in cigar:
        if cigar_op in (0, 7, 8):
            total_matches += num

    return total_matches


class FilteredRecord:
    def __init__(self, args, record, statistics):
        # Flip bit so that 1 represents improper segments
        flag = record.flag ^ BAM_PROPER_SEGMENTS

        insert_size = None
        if flag & BAM_SEGMENTED:
            passes_filters = not flag & args.pe_filter_mask
            insert_size = abs(record.template_length)
        else:
            passes_filters = not flag & args.se_filter_mask

        # Whole flags are counted here and then bulk split into individual types later
        statistics["flags"][flag] += 1

        if record.mapping_quality < args.min_mapping_quality:
            statistics["filters"]["low_mapping_quality"] += 1
            passes_filters = False

        cigartuples = record.cigartuples
        cigarmatches = None
        if cigartuples is not None:
            cigarmatches = count_mapped_bases(tuple(cigartuples))
            if cigarmatches / record.query_length < args.min_mapped_fraction:
                statistics["filters"]["low_mapped_fraction"] += 1
                passes_filters = False

        # Insert size must be reasonable, but is probably not reliable for flagged reads
        if insert_size and passes_filters:
            # Note that insert sizes get counted twice for paired reads here
            json_insert_size = min(args.json_insert_size_cap, insert_size)
            if args.min_insert_size <= insert_size <= args.max_insert_size:
                statistics["insert_sizes"]["passed"][json_insert_size] += 1
            else:
                statistics["insert_sizes"]["failed"][json_insert_size] += 1
                statistics["filters"]["bad_insert_size"] += 1
                passes_filters = False

        self.record = record
        self.statistics = statistics
        self.passes_filters = passes_filters
        self.cigarmatches = cigarmatches

    def finalize(self):
        record = self.record

        key = "passed" if self else "failed"
        self.statistics["totals"][key]["reads"] += 1
        self.statistics["query_lengths"][key][record.query_length] += 1
        self.statistics["mapping_quality"][key][record.mapping_quality] += 1

        if self.cigarmatches is not None:
            self.statistics["matches"][key][self.cigarmatches] += 1
            self.statistics["matches_pct"][key][
                int(self.cigarmatches * 100 / record.query_length)
            ] += 1

        return bool(self), record

    def set_orphan(self):
        self.passes_filters = False
        self.statistics["filters"]["orphan_reads"] += 1

    @property
    def is_mate_1(self):
        flag = self.record.flag
        return (
            (flag & BAM_SEGMENTED)
            and (flag & BAM_IS_FIRST_SEGMENT)
            and not (flag & BAM_IS_LAST_SEGMENT)
        )

    @property
    def is_mate_2(self):
        flag = self.record.flag
        return (
            (flag & BAM_SEGMENTED)
            and (flag & BAM_IS_LAST_SEGMENT)
            and not (flag & BAM_IS_FIRST_SEGMENT)
        )

    def __bool__(self):
        return self.passes_filters


def evaluate_alignments(args, records, statistics):
    # Group statistics by type (merged, paired, etc.)
    statistics = statistics[identify_read(records[0])]

    mate_1_record = None
    mate_2_record = None
    filtered_records = []
    for record in records:
        filtered_record = FilteredRecord(args, record, statistics)

        if filtered_record:
            # FIXME: Check for multiple reads flagged as mate 1 or 2?
            # FIXME: Handle internal segments (flagged as mate 1 AND mate 2)?
            if filtered_record.is_mate_1:
                assert mate_1_record is None
                mate_1_record = filtered_record
            elif filtered_record.is_mate_2:
                assert mate_2_record is None
                mate_2_record = filtered_record

        filtered_records.append(filtered_record)

    if not args.allow_orphan_mates:
        if mate_1_record and not mate_2_record:
            mate_1_record.set_orphan()
        elif not mate_1_record and mate_2_record:
            mate_2_record.set_orphan()

    for record in filtered_records:
        yield record.finalize()


def process_reads_by_queryname(args, handle, statistics):
    for _, records in groupby(BAMTimer(handle), lambda it: it.query_name):
        yield from evaluate_alignments(args, tuple(records), statistics)


def process_individual_reads(args, handle, statistics):
    for record in BAMTimer(handle):
        filtered_record = FilteredRecord(
            args, record, statistics[identify_read(record)]
        )

        yield filtered_record.finalize()


def calculate_flag_statistics(args, statistics):
    for metrics in statistics.values():
        filter_counts = metrics["filters"]

        for flag, count in metrics.pop("flags").items():
            if flag & BAM_SEGMENTED:
                masks = args.named_pe_filters
            else:
                masks = args.named_se_filters

            for label, mask in masks.items():
                if flag & mask:
                    filter_counts[label] += count


def calculate_coverage_statistics(statistics, genome_size):
    for metrics in statistics.values():
        for key, counts in metrics["matches"].items():
            total_matches = 0
            for matches, count in counts.items():
                total_matches += matches * count

            metrics["totals"][key]["bases"] = total_matches
            metrics["totals"][key]["coverage"] = total_matches / genome_size


def normalize_pe_insert_sizes(statistics):
    for group in statistics.values():
        insert_sizes = group["insert_sizes"]
        for key, distribution in insert_sizes.items():
            # Both mates are counted above, so divide by two
            insert_sizes[key] = {
                insert_size: count // 2 for insert_size, count in distribution.items()
            }


def calculate_totals(src, dst):
    for key, value in src.items():
        if isinstance(value, (int, float)):
            dst[key] = dst.get(key, 0) + value
        elif isinstance(value, dict):
            calculate_totals(value, dst.setdefault(key, {}))
        else:
            raise ValueError(key)


def calculate_statistics(args, handle, statistics):
    genome_size = sum(handle.lengths)

    normalize_pe_insert_sizes(statistics)
    calculate_flag_statistics(args, statistics)
    calculate_coverage_statistics(
        statistics=statistics,
        genome_size=genome_size,
    )

    for group, metrics in tuple(statistics.items()):
        totals = metrics["totals"]
        if not any(total["reads"] for total in totals.values()):
            statistics.pop(group)

    totals = {}
    for counts in tuple(statistics.values()):
        calculate_totals(counts, totals)

    statistics["*"] = totals

    return statistics


def write_json(args, in_bam, statistics):
    lengths = in_bam.lengths
    statistics = calculate_statistics(args, in_bam, statistics)

    with args.out_json.open("wt") as handle:
        json.dump(
            {
                "input": str(args.in_bam.absolute()),
                "output_passed": str(args.out_passed.absolute())
                if args.out_passed
                else None,
                "output_failed": str(args.out_failed.absolute())
                if args.out_failed
                else None,
                "statistics": statistics,
                "genome": {
                    "ncontigs": len(lengths),
                    "size": sum(lengths),
                },
            },
            handle,
            sort_keys=True,
        )


def initialize_statistics(args, handle):
    def _passed_and_failed():
        return {
            "passed": defaultdict(int),
            "failed": defaultdict(int),
        }

    readgroups = {}
    for readgroup in handle.header["RG"]:
        readgroups[readgroup["ID"]] = readgroup["DS"]

    # Counts of flags per read group
    statistics = {}
    for group in set(readgroups.values()):
        statistics[group] = {
            "totals": _passed_and_failed(),
            "filters": dict.fromkeys(args.named_pe_filters, 0),
            "query_lengths": _passed_and_failed(),
            "insert_sizes": _passed_and_failed(),
            "matches": _passed_and_failed(),
            "matches_pct": _passed_and_failed(),
            "mapping_quality": _passed_and_failed(),
            "flags": defaultdict(int),
        }

        if args.min_insert_size > 0 or args.max_insert_size < float("inf"):
            statistics[group]["filters"]["bad_insert_size"] = 0

        if args.min_matches > 0:
            statistics[group]["filters"]["few_mapped"] = 0

        if args.min_mapped_fraction > 0:
            statistics[group]["filters"]["low_mapped_fraction"] = 0

        if args.min_mapping_quality > 0:
            statistics[group]["filters"]["low_mapping_quality"] = 0

        if not args.allow_orphan_mates:
            statistics[group]["filters"]["orphan_reads"] = 0

    return statistics, {key: statistics[group] for key, group in readgroups.items()}


def configure_flag_filters(args):
    args.named_se_filters = {
        "unmapped": BAM_READ_IS_UNMAPPED,
        "secondary": BAM_SECONDARY_ALIGNMENT,
        "qc_failed": BAM_QUALITY_CONTROL_FAILED,
        "pcr_duplicate": BAM_PCR_DUPLICATE,
        "supplementary": BAM_SUPPLEMENTARY_ALIGNMENT,
    }

    args.named_pe_filters = dict(args.named_se_filters)

    if not args.allow_orphan_mates:
        args.named_pe_filters["unmapped_mate"] = BAM_NEXT_IS_UNMAPPED

    if not args.allow_improper_pairs:
        # Proper segments is a good thing, so we flip the bit when testing it
        args.named_pe_filters["improper_pair"] = BAM_PROPER_SEGMENTS

    args.se_filter_mask = sum(args.named_se_filters.values())
    args.pe_filter_mask = sum(args.named_pe_filters.values())


def is_queryname_ordering_required(args):
    return (
        args.min_insert_size > 0
        or args.max_insert_size < float("inf")
        or args.min_matches > 0
        or args.min_mapped_fraction > 0
        or args.min_mapping_quality > 0
    ) and not args.allow_orphan_mates


def is_queryname_ordered(handle):
    return handle.header.get("HD", {}).get("SO") == "queryname"


def parse_args(argv):
    parser = ArgumentParser("paleomix ngs:finalize_bam")
    parser.add_argument("in_bam", type=Path, default=Path("-"), nargs="?")

    parser.add_argument(
        "--threads",
        default=4,
        type=int,
        help="Number of threads used for reading-writing BAM files",
    )

    group = parser.add_argument_group("Output")
    group.add_argument(
        "--out-passed",
        type=Path,
        help="Output path for BAM reads that pass all filters",
    )
    group.add_argument(
        "--out-passed-uncompressed",
        action="store_true",
        help="Output uncompressed BAM for reads that pass all filters",
    )
    group.add_argument(
        "--out-failed",
        type=Path,
        help="Output path for BAM reads that fail one or more filters",
    )
    group.add_argument(
        "--out-json",
        type=Path,
        help="Output path for filtering/alignment statistics in JSON format",
    )

    group = parser.add_argument_group("Filters")
    group.add_argument(
        "--min-insert-size",
        type=int,
        default=0,
        help="Minimum insert size for paired/merged records",
    )
    group.add_argument(
        "--max-insert-size",
        type=int,
        default=float("inf"),
        help="Maximum insert size for paired/merged reads",
    )
    group.add_argument(
        "--min-matches",
        type=int,
        default=0,
        help="Minimum number of aligned bases (cigar M, =, and X)",
    )
    group.add_argument(
        "--min-mapping-quality",
        type=int,
        default=0,
        help="Minimum mapping quality for alignments",
    )
    group.add_argument(
        "--min-mapped-fraction",
        type=float,
        default=0,
        help="Minimum fraction (0.0-1.0) of aligned bases (M, =, X) in an alignment",
    )

    group.add_argument(
        "--allow-improper-pairs",
        action="store_true",
        help="If set, the proper-pairs bit (0x2) is not required for paired reads",
    )
    group.add_argument(
        "--allow-orphan-mates",
        action="store_true",
        help="If set, one read in a pair being filtered does not cause the other read "
        "to be filtered",
    )

    group = parser.add_argument_group("JSON")
    group.add_argument(
        "--json-insert-size-cap",
        type=int,
        default=1_001,
        help="Cap insert sizes at this value for reporting purposes. "
        "Set to <= 0 to disable. Does not affect filtering",
    )
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    initialize_console_logging()

    if args.json_insert_size_cap <= 0:
        args.json_insert_size_cap = float("inf")

    log = logging.getLogger("ngs:finalize_bam")
    if not (args.out_passed or args.out_failed or args.out_json):
        log.error("No --out-* arguments; please specify at least one output file")
        return 1

    # Decide on which filters to enable/count in statistics
    configure_flag_filters(args)

    log.info("Reading alignents from %s", args.in_bam)
    in_bam = pysam.AlignmentFile(str(args.in_bam), threads=args.threads)

    # Some filters require that we consider both the mate 1 and the mate 2 read at the
    # same time, in which case the file has to be query-name ordered. This is because
    # some filters may apply to one read but not the other (e.g. mapping quality).
    if is_queryname_ordered(in_bam):
        bamfile_reader = process_reads_by_queryname
    elif is_queryname_ordering_required(args):
        log.error("A queryname sorted input BAM is required for these filters!")
        return 1
    else:
        bamfile_reader = process_individual_reads

    if args.out_passed:
        log.info("Writing proper alignents to %s", args.out_passed)
        out_passed = pysam.AlignmentFile(
            args.out_passed, "wb", template=in_bam, threads=args.threads
        )
    else:
        out_passed = None

    if args.out_failed:
        log.info("Writing failed alignents/unaligned reads to %s", args.out_failed)
        out_failed = pysam.AlignmentFile(
            args.out_failed, "wb", template=in_bam, threads=args.threads
        )
    else:
        out_failed = None

    statistics, readgroup_to_statistics = initialize_statistics(args, in_bam)
    for passed, record in bamfile_reader(args, in_bam, readgroup_to_statistics):
        if passed:
            if out_passed:
                out_passed.write(record)
        else:
            if out_failed:
                out_failed.write(record)

    if out_passed:
        out_passed.close()

    if out_failed:
        out_failed.close()

    if args.out_json:
        write_json(args, in_bam, statistics)

    in_bam.close()

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv[1:]))
