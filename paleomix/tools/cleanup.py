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
import sys

import pysam

import paleomix.tools.factory

import paleomix.common.argparse as argparse
import paleomix.common.procs as processes


# Mask to select flags that are relevant to SE reads; this excludes flags where
# no assumptions can if 0x1 is not set, per the SAM specification (see below).
_SE_FLAGS_MASK = ~(0x2 | 0x8 | 0x20 | 0x40 | 0x80)


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


def _cleanup_record(record):
    """Cleans up the properties of a BAM record, ensuring that only appropriate
    flags and values are set, such that the record follows section 1.4 of the
    SAM specification (https://samtools.github.io/hts-specs/SAMv1.pdf). The
    record itself (or a new record) is returned.
    """
    if not record.is_paired:
        # Unset 0x2 (properly aligned), 0x8 (next mate unmapped),
        # 0x20 (next mate reverse), 0x40 (first mate), and 0x80 (last mate).
        record.flag = record.flag & (~0xEA)
        record.rnext = -1
        record.pnext = -1
        record.tlen = 0
    elif record.mate_is_unmapped and record.has_tag("MC"):
        # Picard ValidateSamFile (2.9.1) objects to MC tags for unmapped mates,
        # which are currently added by SAMTools (v1.4).
        tags = record.get_tags(with_value_type=True)
        record.set_tags([tag for tag in tags if tag[0] != "MC"])

    if record.is_unmapped:
        record.mapq = 0
        record.cigar = None
        # Unset 0x2 (properly aligned), 0x100 (secondary), and 0x800 (chimeric)
        record.flag = record.flag & (~0x902)

        if record.mate_is_unmapped:
            record.rnext = -1
            record.pnext = -1

        # Per the spec, unmapped reads should be placed with their mate
        record.tid = record.rnext
        record.pos = record.pnext
        record.tlen = 0
    elif record.mate_is_unmapped:
        record.rnext = record.tid
        record.pnext = record.pos
        record.tlen = 0

    return record


def _filter_record(args, record):
    """Returns True if the record should be filtered (excluded), based on the
    --exclude-flags and --require-flags options. Certain flags are ignored when
    filtering SE reads, namely those not included in _SE_FLAGS_MASK (above).
    """
    if record.is_paired:
        exclude_flags = args.exclude_flags
        require_flags = args.require_flags
    else:
        exclude_flags = args.exclude_flags & _SE_FLAGS_MASK
        require_flags = args.require_flags & _SE_FLAGS_MASK

    if record.flag & exclude_flags:
        return True
    elif ~(record.flag & require_flags) & require_flags:
        return True

    return False


def _cleanup_unmapped(args):
    """Reads a BAM (or SAM, if cleanup_sam is True) file from STDIN, and
    filters reads according to the filters specified in the commandline
    arguments 'args'. The resulting records are written to STDOUT in
    uncompressed BAM format. The output BAM is marked as sorted (under the
    assumption that 'samtools sort' is to be run on the output) and PG tags are
    updated if specified in the args.
    """

    filter_by_flag = bool(args.exclude_flags or args.require_flags)
    with pysam.AlignmentFile("-") as input_handle:
        header = dict(input_handle.header)
        _set_sort_order(header)
        _set_pg_tags(header, args.update_pg_tag)
        if args.rg_id is not None:
            _set_rg_tags(header, args.rg_id, args.rg)

        with pysam.AlignmentFile("-", "wbu", header=header) as output_handle:
            for record in input_handle:
                # Ensure that the properties make sense before filtering
                record = _cleanup_record(record)

                if not record.is_unmapped and (record.mapq < args.min_quality):
                    continue
                elif filter_by_flag and _filter_record(args, record):
                    continue

                if args.rg_id is not None:
                    # Ensure that only one RG tag is set
                    tags = record.get_tags(with_value_type=True)
                    tags = [tag for tag in tags if tag[0] != "RG"]
                    tags.append(("RG", args.rg_id, "Z"))
                    record.set_tags(tags)

                output_handle.write(record)

    return 0


def _build_wrapper_command(args):
    bam_cleanup = paleomix.tools.factory.new("cleanup")
    if args.fasta is not None:
        bam_cleanup.set_option("--fasta", args.fasta)
    bam_cleanup.set_option("--temp-prefix", args.temp_prefix)
    bam_cleanup.set_option("--min-quality", str(args.min_quality))
    bam_cleanup.set_option("--exclude-flags", hex(args.exclude_flags))

    for value in args.update_pg_tag:
        bam_cleanup.add_option("--update-pg-tag", value)

    if args.rg_id is not None:
        bam_cleanup.set_option("--rg-id", args.rg_id)
        for value in args.rg:
            bam_cleanup.add_option("--rg", value)

    return bam_cleanup.call


def _run_cleanup_pipeline(args):
    bam_cleanup = _build_wrapper_command(args)
    commands = []

    try:
        if args.paired_end:
            # Convert input to (uncompressed) BAM and fix mate information for PE reads
            commands.append(["samtools", "fixmate", "-O", "bam", "-", "-"])

        # Cleanup / filter reads. Must be done after 'fixmate', as BWA may produce
        # hits where the mate-unmapped flag is incorrect, which 'fixmate' fixes.
        commands.append(bam_cleanup + ["cleanup"])

        # Sort by coordinates and output uncompressed BAM
        commands.append(
            ["samtools", "sort", "-l", "0", "-O", "bam", "-T", args.temp_prefix]
        )

        # Update NM and MD tags; output BAM (-b) to stdout
        if args.fasta is not None:
            commands.append(["samtools", "calmd", "-b", "-", args.fasta])

        procs = []
        last_out = sys.stdin
        for cmd in commands:
            proc_stdout = None if cmd is commands[-1] else processes.PIPE
            procs.append(processes.open_proc(cmd, stdin=last_out, stdout=proc_stdout))

            last_out.close()
            last_out = procs[-1].stdout

        return int(any(processes.join_procs(procs)))
    except Exception:
        for proc in procs.values():
            proc.terminate()
        raise


def parse_args(argv):
    """Parses a list of command-line arguments, excluding the program name."""
    prog = "paleomix cleanup"
    usage = "%s --temp-prefix prefix --fasta reference.fasta < in.sam" % (prog,)

    parser = argparse.ArgumentParser(prog=prog, usage=usage)
    # "Hidden" commands, invoking the various sub-parts of this script
    parser.add_argument("command", nargs="?", help=argparse.SUPPRESS)

    parser.add_argument(
        "--fasta",
        default=None,
        help="Reference FASTA sequence; if set, the calmd command is used to "
        "re-calculate MD tags.",
    )
    parser.add_argument(
        "--temp-prefix", required=True, help="REQUIRED: Prefix for temp files"
    )
    parser.add_argument(
        "-q",
        "--min-quality",
        type=int,
        default=0,
        help="Exclude aligned reads with a mapping quality below this value; note "
        "that this filter ONLY applies to aligned reads",
    )
    parser.add_argument(
        "-f",
        "--require-flags",
        default=0,
        type=lambda value: int(value, 0),  # Handle hex, etc.
        help="Only include reads with all of these flags set; note that flags only "
        "valid for paired-end reads (0x2, 0x8, 0x20, 0x40, 0x80) are ignored when "
        "processing single-end reads.",
    )
    parser.add_argument(
        "-F",
        "--exclude-flags",
        default=0,
        type=lambda value: int(value, 0),  # Handle hex, etc.
        help="Exclude reads with any of these flags set; note that flags only valid "
        "for paired-end reads (0x2, 0x8, 0x20, 0x40, 0x80) are ignored when "
        "processing single-end reads.",
    )
    parser.add_argument(
        "--paired-end",
        default=False,
        action="store_true",
        help="If enabled, additional processing of PE reads is carried out, including "
        "updating of mate information [Default: off]",
    )

    parser.add_argument(
        "--update-pg-tag",
        default=[],
        action="append",
        help="Update one PG tags with the given values, creating the tag if it does "
        'not already exist. Takes arguments in the form "PGID:TAG:VALUE".',
    )
    parser.add_argument(
        "--rg-id",
        default=None,
        help="If set, the read-group is overwritten based on tags set using the --rg "
        "option, using the id specified using --rg-id.",
    )
    parser.add_argument(
        "--rg",
        default=[],
        action="append",
        help="Readgroup in the form 'ID:TAG:VALUE'.",
    )

    args = parser.parse_args(argv)
    if args.command not in (None, "cleanup"):
        parser.error("unrecognized arguments: %s" % (args.command,))

    return args


def main(argv):
    """Main function; returns 0 on success, non-zero otherwise."""
    args = parse_args(argv)

    if args.command == "cleanup":
        return _cleanup_unmapped(args)
    elif args.command:
        raise NotImplementedError("Unexpected command %r" % (args.command,))

    sys.stderr.write("Reading SAM file from STDIN\n")
    return _run_cleanup_pipeline(args)
