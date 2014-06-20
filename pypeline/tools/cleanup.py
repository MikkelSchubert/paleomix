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
"""
Reads SAM or BAM from STDIN and outputs cleaned and sorted BAM to STDOUT.

The cleanup involves setting header-flags, fixing mate information, clearing
various fields of unmapped reads (due to BWA leaving CIGARs in unmapped reads),
sorting the input, and updating MD and NM tags.

This script also solves a problem with parsing record-less files by SAMTools,
which fails with a parse error if the input SAM file only contains a header.
This is mainly a problem when aligning ancient DNA lanes against mt genomes,
as there may not be any hits in an entire lane. For example, the following
command will not work:
$ samtools view -H INPUT.BAM | samtools view -Sbu -

"""
import os
import sys
import copy
import time
import argparse
import subprocess

import pysam


def _set_sort_order(header):
    """Updates a BAM header to indicate coordinate sorting."""
    hd_dict = header.setdefault("HD", {"GO": "none", "VN": "1.0"})
    hd_dict["SO"] = "coordinate"


def _set_pg_tags(header, tags):
    """Updates PG tags in a BAM header, taking a sequence of ID:TAG:VALUEs."""
    for tag in tags:
        pg_id, pg_field, pg_value = tag.split(":")

        for pg_dict in header.setdefault("PG", []):
            if pg_dict.get("ID") == pg_id:
                pg_dict[pg_field] = pg_value
                break
        else:
            header["PG"].append({"ID": pg_id, pg_field: pg_value})


def _set_rg_tags(header, rg_id, rg_tags):
    """Updates RG tags in a BAM header, taking a sequence of TAG:VALUEs."""
    readgroup = {"ID": rg_id}
    for tag in rg_tags:
        rg_field, rg_value = tag.split(":")
        readgroup[rg_field] = rg_value
    header["RG"] = [readgroup]



def _pipe_to_bam():
    """Simply pipes a BAM/SAM file to stdout; this is required to handle SAM
    files that do not contain records (i.e. only a header), which are not
    properly handled by "samtools view -S -", resulting in a parse failure.
    """
    with pysam.Samfile("-", "r") as input_handle:
        with pysam.Samfile("-", "wbu", template=input_handle) as output_handle:
            for record in input_handle:
                output_handle.write(record)

    return 0


def _cleanup_record(record):
    """Marks a BAM record as unmapped, clearing relevant fields and/or setting
    fields to match those of the mate (if mapped). An updated (possibly new)
    record is returned.
    """
    if record.cigar:
        # Build a new read; this is nessesary, as it is not possible
        # to clean the CIGAR string on an existing record in current
        # versions of Pysam.
        unmapped_read = pysam.AlignedRead()
        unmapped_read.qname = record.qname
        unmapped_read.flag = record.flag
        unmapped_read.seq = record.seq
        unmapped_read.qual = record.qual
        unmapped_read.tags = record.tags

        if not record.mate_is_unmapped:
            unmapped_read.rnext = record.rnext
            unmapped_read.pnext = record.pnext
        else:
            unmapped_read.rnext = -1
            unmapped_read.pnext = -1
        unmapped_read.tid = unmapped_read.rnext

        # Set .pos TWICE; this is a workaround for a bug in current versions
        # of pysam, in which the bin in the record is re-calculated BEFORE
        # the new position value is set, using the old pos value.
        unmapped_read.pos = unmapped_read.pnext  # Update 1 of 2
        unmapped_read.pos = unmapped_read.pnext  # Update 2 of 2

        return unmapped_read
    else:
        record.mapq = 0
        if record.mate_is_unmapped:
            record.rnext = -1
            record.pnext = -1
        record.tid = record.rnext
        record.pos = record.pnext
        record.tlen = 0

        return record


def _cleanup_unmapped(args, cleanup_sam):
    """Reads a BAM (or SAM, if cleanup_sam is True) file from STDIN, and
    filters reads according to the filters specified in the commandline
    arguments 'args'. The resulting records are written to STDOUT in
    uncompressed BAM format. The output BAM is marked as sorted (under the
    assumption that 'samtools sort' is to be run on the output) and PG tags are
    updated if specified in the args.
    """
    spec = "r" if cleanup_sam else "rb"
    with pysam.Samfile("-", spec) as input_handle:
        header = copy.deepcopy(input_handle.header)
        _set_sort_order(header)
        _set_pg_tags(header, args.update_pg_tag)
        if args.rg_id is not None:
            _set_rg_tags(header, args.rg_id, args.rg)

        with pysam.Samfile("-", "wbu", header=header) as output_handle:
            for record in input_handle:
                if (record.mapq < args.min_quality) \
                        or (record.flag & args.exclude_flags):
                    continue

                if record.is_unmapped:
                    # Unmapped read; clear all non-required fields
                    record = _cleanup_record(record)
                elif record.mate_is_unmapped:
                    # Unmapped mate
                    record.rnext = record.tid
                    record.pnext = record.pos
                    record.tlen = 0

                if args.rg_id is not None:
                    # Ensure that only one RG tag is set
                    tags = [(key, value) for (key, value) in record.tags
                            if key != "RG"]
                    tags.append(("RG", args.rg_id))
                    record.tags = tags

                output_handle.write(record)

    return 0


def _join_procs(procs):
    """Joins a set of Popen processes. If a processes fail, the remaining
    processes are terminated. The function returns 1 if any processes failed,
    0 otherwise.
    """
    sleep_time = 0.05
    commands = list(enumerate(procs))
    return_codes = [None] * len(commands)
    while commands:
        for (index, command) in list(commands):
            if command.poll() is not None:
                return_codes[index] = command.wait()
                commands.remove((index, command))
                sleep_time = 0.05
            elif any(return_codes):
                command.terminate()
                return_codes[index] = command.wait()
                commands.remove((index, command))
                sleep_time = 0.05

        time.sleep(sleep_time)
        sleep_time = min(1, sleep_time * 2)

    return 1 if any(return_codes) else 0


def _setup_single_ended_pipeline(procs, bam_cleanup):
    # Convert input to BAM and cleanup / filter reads
    call_pipe = bam_cleanup + ['cleanup-sam']

    # Convert input to (uncompressed) BAM
    procs["pipe"] = subprocess.Popen(call_pipe,
                                     stdin=sys.stdin,
                                     stdout=subprocess.PIPE,
                                     close_fds=True)
    sys.stdin.close()

    return procs["pipe"]


def _setup_paired_ended_pipeline(procs, bam_cleanup):
    call_pipe = bam_cleanup + ['pipe']

    # Convert input to (uncompressed) BAM
    procs["pipe"] = subprocess.Popen(call_pipe,
                                     stdin=sys.stdin,
                                     stdout=subprocess.PIPE,
                                     close_fds=True)
    sys.stdin.close()

    # Fix mate information for PE reads
    call_fixmate = ['samtools', 'fixmate', '-', '-']
    procs["fixmate"] = subprocess.Popen(call_fixmate,
                                        stdin=procs["pipe"].stdout,
                                        stdout=subprocess.PIPE,
                                        close_fds=True)
    procs["pipe"].stdout.close()

    # Cleanup / filter reads. Must be done after 'fixmate', as BWA may produce
    # hits where the mate-unmapped flag is incorrect, which 'fixmate' fixes.
    procs["cleanup"] = subprocess.Popen(bam_cleanup + ['cleanup'],
                                        stdin=procs["fixmate"].stdout,
                                        stdout=subprocess.PIPE,
                                        close_fds=True)
    procs["fixmate"].stdout.close()

    return procs["cleanup"]


def _run_cleanup_pipeline(args):
    # Ensure that the same script is used, even if it is not in the PATH
    bam_cleanup = [sys.argv[0]]
    if os.path.basename(bam_cleanup[0]) == "paleomix":
        bam_cleanup.append("cleanup")

    bam_cleanup.extend(('--fasta', args.fasta,
                        '--temp-prefix', args.temp_prefix,
                        '--min-quality', str(args.min_quality),
                        '--exclude-flags', hex(args.exclude_flags)))

    for value in args.update_pg_tag:
        bam_cleanup.extend(('--update-pg-tag', value))

    if args.rg_id is not None:
        bam_cleanup.extend(('--rg-id', args.rg_id))
        for value in args.rg:
            bam_cleanup.extend(('--rg', value))

    procs = {}
    try:
        # Update 'procs' and get the last process in the pipeline
        if args.paired_ended:
            last_proc = _setup_paired_ended_pipeline(procs, bam_cleanup)
        else:
            last_proc = _setup_single_ended_pipeline(procs, bam_cleanup)

        # Sort, output to stdout (-o)
        call_sort = ['samtools', 'sort', "-o", "-", args.temp_prefix]
        procs["sort"] = subprocess.Popen(call_sort,
                                         stdin=last_proc.stdout,
                                         stdout=subprocess.PIPE,
                                         close_fds=True)
        last_proc.stdout.close()

        # Update NM and MD tags; output BAM (-b) to stdout
        call_calmd = ['samtools', 'calmd', '-b', '-', args.fasta]
        procs["calmd"] = subprocess.Popen(call_calmd,
                                          stdin=procs["sort"].stdout)
        procs["sort"].stdout.close()

        return _join_procs(procs.values())
    except:
        for proc in procs.itervalues():
            proc.terminate()
        raise


def parse_args(argv):
    parser = argparse.ArgumentParser()
    # "Hidden" commands, invoking the various sub-parts of this script
    parser.add_argument('command', choices=('pipe', 'cleanup', 'cleanup-sam'),
                        nargs="?", help=argparse.SUPPRESS)
    # Specifies if the 'cleanup' step should expect SAM input
    parser.add_argument('--cleanup-sam', default=False, action="store_true",
                        help=argparse.SUPPRESS)

    parser.add_argument('--fasta', help="REQUIRED: Reference FASTA sequence",
                        required=True)
    parser.add_argument('--temp-prefix', required=True,
                        help="REQUIRED: Prefix for temp files")
    parser.add_argument("-q", "--min-quality", type=int, default=0,
                        help="Equivalent to \"samtools view -q\"  "
                             "\n[Default: %(default)s]")
    parser.add_argument("-F", "--exclude-flags", default=0,
                        type=lambda value: int(value, 0),  # Handle hex, etc.
                        help="Equivalent to \"samtools view -F\" "
                             "[Default: %(default)s]")
    parser.add_argument('--paired-ended', default=False, action="store_true",
                        help='If enabled, additional processing of PE reads'
                             'is carried out, including updating of mate '
                             'information [Default: off]')

    parser.add_argument("--update-pg-tag", default=[], action="append",
                        help="Update one PG tags with the given values, "
                             "creating the tag if it does not already exist. "
                             "Takes arguments in the form \"PGID:TAG:VALUE\".")
    parser.add_argument('--rg-id', default=None,
                        help="If set, the read-group is overwritten based "
                             "on tags set using the --rg option, using the "
                             "id specified using --rg-id.")
    parser.add_argument('--rg', default=[], action="append",
                        help="Create readgroup with ") # TODO

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    if args.command == "pipe":
        return _pipe_to_bam()
    elif args.command == "cleanup":
        return _cleanup_unmapped(args, cleanup_sam=False)
    elif args.command == "cleanup-sam":
        return _cleanup_unmapped(args, cleanup_sam=True)

    sys.stderr.write("Reading SAM file from STDIN ...\n")
    return _run_cleanup_pipeline(args)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
