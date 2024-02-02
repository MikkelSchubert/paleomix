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
import logging
import os
from dataclasses import dataclass, field
from typing import Callable

import pysam

from paleomix.common import argparse
from paleomix.common.fileutils import swap_ext
from paleomix.common.formats.bed import BEDRecord, read_bed_file, sort_bed_by_bamfile


class BAMStatsError(RuntimeError):
    pass


@dataclass
class Args:
    infile: str
    outfile: str
    target_name: str
    regions_fpath: str | None
    max_contigs: int
    ignore_readgroups: bool
    overwrite_output: bool
    get_readgroup_func: Callable[[pysam.AlignedSegment], str | None]
    regions: list[BEDRecord] = field(default_factory=list)


@dataclass
class Readgroup:
    key: str | None
    sample: str
    library: str


def collect_readgroups(
    args: Args,
    handle: pysam.AlignmentFile,
) -> list[Readgroup]:
    readgroups = [Readgroup(key=None, sample="<NA>", library="<NA>")]
    if args.ignore_readgroups:
        return readgroups

    for readgroup in handle.header.get("RG", ()):
        readgroups.append(
            Readgroup(
                key=readgroup["ID"],
                sample=readgroup["SM"],
                library=readgroup["LB"],
            )
        )

    return readgroups


def collect_references(args: Args, handle: pysam.AlignmentFile) -> dict[str, int]:
    if args.regions:
        lengths: dict[str, int] = collections.defaultdict(int)
        for region in args.regions:
            lengths[region.name] += region.end - region.start

        lengths = dict(lengths)
    elif handle.nreferences <= args.max_contigs:
        lengths = dict(zip(handle.references, handle.lengths))
    else:
        lengths = {"<Genome>": sum(handle.lengths)}

    return lengths


def collect_bed_regions(filename: str) -> list[BEDRecord]:
    cache: dict[str, str] = {}
    regions: list[BEDRecord] = []
    for record in read_bed_file(filename):
        name = record.name
        if not name:
            name = f"{record.contig}*"

        record.contig = cache.setdefault(record.contig, record.contig)
        record.name = cache.setdefault(name, name)

        regions.append(
            BEDRecord(
                contig=record.contig,
                start=record.start,
                end=record.end,
                name=name,
            )
        )

    return regions


def parse_arguments(argv: list[str], ext: str) -> Args:
    prog = "paleomix {}".format(ext.strip("."))
    usage = f"{prog} [options] sorted.bam [out{ext}]"
    parser = argparse.ArgumentParser(prog=prog, usage=usage)

    parser.add_argument(
        "infile",
        metavar="BAM",
        help="Filename of a sorted BAM file. If set to '-' "
        "the file is read from STDIN.",
    )
    parser.add_argument(
        "outfile",
        metavar="OUTPUT",
        nargs="?",
        help="Filename of output table; defaults to name of "
        f"the input BAM with a '{ext}' extension. If "
        "set to '-' the table is printed to STDOUT.",
    )
    parser.add_argument(
        "--target-name",
        default=None,
        metavar="NAME",
        help="Name used for 'Target' column; defaults to the "
        "filename of the BAM file.",
    )
    parser.add_argument(
        "--regions-file",
        default=None,
        dest="regions_fpath",
        help="BED file containing regions of interest; {} "
        "is calculated only for these grouping by the "
        "name used in the BED file, or the contig name "
        "if no name has been specified for a record.".format(ext.strip(".")),
    )
    parser.add_argument(
        "--max-contigs",
        default=100,
        type=int,
        help="The maximum number of contigs allowed in a BAM "
        "file. If this number is exceeded, the entire "
        "set of contigs is aggregated into one pseudo-"
        "contig named '<Genome>'. This is done to "
        "limit table sizes",
    )
    parser.add_argument(
        "--ignore-readgroups",
        default=False,
        action="store_true",
        help="Ignore readgroup information in reads, and only "
        "provide aggregated statistics; this is required "
        "if readgroup information is missing or partial",
    )
    parser.add_argument(
        "--overwrite-output",
        default=False,
        action="store_true",
        help="Overwrite output file if it it exists; by "
        "default, the script will terminate if the file "
        "already exists.",
    )

    args = parser.parse_args(argv)
    if not args.outfile:
        args.outfile = swap_ext(args.infile, ext)

    if args.ignore_readgroups:
        get_readgroup_func = _get_readgroup_ignored
    else:
        get_readgroup_func = _get_readgroup

    if not args.target_name:
        if args.infile == "-":
            args.target_name = "<STDIN>"
        else:
            args.target_name = os.path.basename(args.infile)

    if os.path.exists(args.outfile) and not args.overwrite_output:
        parser.error(
            f"Destination filename already exists ({args.outfile!r}); use option "
            "--overwrite-output to allow overwriting of this file."
        )

    return Args(
        infile=args.infile,
        outfile=args.outfile,
        target_name=args.target_name,
        regions_fpath=args.regions_fpath,
        max_contigs=args.max_contigs,
        ignore_readgroups=args.ignore_readgroups,
        overwrite_output=args.overwrite_output,
        get_readgroup_func=get_readgroup_func,
    )


def main_wrapper(
    process_func: Callable[[Args, pysam.AlignmentFile], int],
    argv: list[str],
    ext: str,
) -> int:
    log = logging.getLogger(__name__)
    args = parse_arguments(argv, ext)
    if args.regions_fpath:
        try:
            args.regions = collect_bed_regions(args.regions_fpath)
        except ValueError as error:
            log.error("Failed to parse BED file %r: %s", args.regions_fpath, error)  # noqa: TRY400
            return 1

    log.info("Opening %r", args.infile)
    with pysam.AlignmentFile(args.infile) as handle:
        header = handle.header.to_dict()
        sort_order = header.get("HD", {}).get("SO")
        if sort_order is None:
            log.warning("BAM file %r is not marked as sorted!", args.infile)
        elif sort_order != "coordinate":
            log.error(
                "BAM file %r is %s-sorted, but coordinate-sorting is required",
                args.infile,
                sort_order,
            )
            return 1

        sort_bed_by_bamfile(handle, args.regions)
        return process_func(args, handle)


def _get_readgroup(record: pysam.AlignedSegment) -> str | None:
    try:
        return record.get_tag("RG")
    except KeyError:
        return None


def _get_readgroup_ignored(_: pysam.AlignedSegment) -> str | None:
    return None
