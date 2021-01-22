#!/usr/bin/env python3
# -*- coding: utf8 -*-
import collections
import heapq
import shlex
import sys

import pysam

from paleomix.common.argparse import ArgumentParser, SUPPRESS
from paleomix.common.sequences import reverse_complement


class PipelineHandler:
    def __init__(self, input_files, output_file, quiet):
        self.duplicate_reads = 0
        self.input_files = input_files
        self.output_file = output_file
        self.quiet = quiet

    def __call__(self, chrom, pos, records, name, seq, qual):
        self.duplicate_reads += 1
        if self.quiet or self.duplicate_reads > 1:
            return

        print("Read was found multiple times:")
        print("    Name:      ", name)
        print("    Sequence:  ", seq)
        print("    Qualities: ", qual)
        print("    Contig:    ", chrom)
        print("    Position:  ", pos)
        print()

        print("Read was found in these BAM files:")
        for filename, records in sorted(records.items()):
            print("   %s in %r" % (_summarize_reads(records), filename))

        print()
        print()

    def finalize(self):
        if not self.duplicate_reads:
            print("No duplicate reads found.")
            return True

        print("A total of", self.duplicate_reads, "reads with duplicates were found.")
        print()
        print("This indicates that the same data has been included multiple times in")
        print("the project. This can be because multiple copies of the same files were")
        print("used, or because one or more files contain multiple copies of the same ")
        print("reads.")
        print()
        print("The 'paleomix dupcheck' command may be used to review the potentially ")
        print("duplicated data in these BAM files:")
        print()
        print("    $ paleomix dupcheck", " ".join(map(shlex.quote, self.input_files)))
        print()
        print()
        print("If this error was a false positive, then you may execute the following ")
        print("command to mark this test as having succeeded:")
        print()
        print("    $ touch", shlex.quote(self.output_file))

        return False


class CLIHandler:
    def __init__(self, quiet=False):
        self.quiet = quiet
        self.duplicate_reads = 0

    def __call__(self, chrom, pos, records, name, seq, qual):
        self.duplicate_reads += 1
        if self.quiet:
            return

        print("%s:%i -- %s %s %s:" % (chrom, pos, name, seq, qual))
        for filename, records in sorted(records.items()):
            print("    - %s:" % (filename,))

            for idx, record in enumerate(records, start=1):
                print("% 8i. " % (idx,), end="")

                if record.is_paired:
                    if record.is_read1:
                        print("Mate 1 read", end="")
                    elif record.is_read2:
                        print("Mate 2 read", end="")
                    else:
                        print("Unpaired read", end="")
                else:
                    print("Unpaired read", end="")

                try:
                    print(" with read-group %r" % (record.get_tag("RG"),), end="")
                except KeyError:
                    pass

                print()

    def finalize(self):
        if self.quiet:
            print("%i" % (self.duplicate_reads,))
        else:
            print("Found %i record(s) with duplicates." % (self.duplicate_reads,))

        return not self.duplicate_reads


def _open_samfiles(handles, filenames):
    sequences = []
    for filename in filenames:
        handle = pysam.AlignmentFile(filename)
        handles.append(handle)

        sequences.append(_read_samfile(handle, filename))

    return heapq.merge(*sequences, key=_key_by_tid_pos)


def _read_samfile(handle, filename):
    for record in handle:
        if record.is_unmapped and (not record.pos or record.mate_is_unmapped):
            # Ignore unmapped reads except when these are sorted
            # according to the mate position (if mapped)
            continue
        elif record.flag & 0x900:
            # Ignore supplementary / secondary alignments
            continue

        yield (record, filename)


def _process_bam_reads(observed_reads, references, position, handler):
    for records_and_filenames in observed_reads.values():
        if len(records_and_filenames) == 1:
            # Most read-names should be obseved at most once at a position
            continue

        result = collections.defaultdict(list)
        for record, filename in records_and_filenames:
            key = (record.is_reverse, record.qname, record.seq, record.qual)
            result[key].append((filename, record))

        for (is_reverse, name, seq, qual), filenames in result.items():
            if len(filenames) == 1:
                # Two reads had same name, but different characterstics
                continue

            records = collections.defaultdict(list)
            for filename, record in filenames:
                records[filename].append(record)

            if is_reverse:
                seq = reverse_complement(seq)
                qual = qual[::-1]

            chrom = references[position[1]]
            pos = position[0]

            handler(chrom, pos, records, name, seq, qual)


def _summarize_reads(records):
    counts = {"mate 1": 0, "mate 2": 0, "unpaired": 0}

    for record in records:
        if record.is_paired:
            if record.is_read1:
                counts["mate 1"] += 1
            elif record.is_read2:
                counts["mate 2"] += 1
            else:
                counts["unpaired"] += 1
        else:
            counts["unpaired"] += 1

    result = []
    for key, value in sorted(counts.items()):
        if value > 1:
            result.append("%i %s reads" % (value, key))
        elif value:
            result.append("%i %s read" % (value, key))

    return ", ".join(result) or "No reads"


def _key_by_tid_pos(record):
    return (record[0].tid, record[0].pos)


def parse_args(argv):
    parser = ArgumentParser(
        prog="paleomix dupcheck",
        description="Attempt to detect reads included multiple times as input based"
        "on the presence of reads with identical names AND sequences. This is "
        "compromise between sensitivity, specificity, and running time.",
    )
    parser.add_argument("files", nargs="+", help="One or more input BAM files.")
    parser.add_argument(
        "--quiet",
        default=False,
        action="store_true",
        help="Only print the number of BAM records where one or more potential "
        "duplicates were identified.",
    )

    # Pipeline mode, showing output/error messages useful for that use-case
    parser.add_argument("--pipeline-output", help=SUPPRESS)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    if args.pipeline_output is not None:
        handler = PipelineHandler(
            quiet=args.quiet,
            input_files=args.files,
            output_file=args.pipeline_output,
        )
    else:
        handler = CLIHandler(quiet=args.quiet)

    handles = []
    last_pos = None
    observed_reads = collections.defaultdict(list)
    reads_iter = _open_samfiles(handles, args.files)
    references = handles[0].references

    for (record, filename) in reads_iter:
        curr_pos = (record.pos, record.tid)
        if curr_pos != last_pos:
            _process_bam_reads(observed_reads, references, last_pos, handler)
            observed_reads.clear()
            last_pos = curr_pos

            # Stop once the trailing, unmapped reads are reached
            if record.tid == -1:
                break

        observed_reads[record.qname].append((record, filename))

    _process_bam_reads(observed_reads, references, last_pos, handler)

    if not handler.finalize():
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
