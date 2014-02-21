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
import os
import sys
import argparse
import collections

import pysam

import pypeline.common.system

from pypeline.common.fileutils import \
    swap_ext
from pypeline.common.bedtools import \
    sort_bed_by_bamfile, \
    read_bed_file


class BAMStatsError(RuntimeError):
    pass


def collect_readgroups(args, handle):
    if args.ignore_readgroups:
        return {None: {"SM": "<NA>", "LB": "<NA>"}}

    readgroups = {}
    for readgroup in handle.header.get("RG", ()):
        key_id = readgroup["ID"]
        sample = readgroup["SM"]
        library = readgroup["LB"]

        readgroups[key_id] = {"SM": sample, "LB": library}
    return readgroups


def collect_references(args, handle):
    if args.regions:
        lengths = collections.defaultdict(int)
        for region in args.regions:
            lengths[region.name] += region.end - region.start

        lengths = dict(lengths)
    elif len(handle.references) <= args.max_contigs:
        lengths = dict(zip(handle.references, handle.lengths))
    else:
        lengths = {"<Genome>": sum(handle.lengths)}

    return lengths


def collect_bed_regions(filename):
    regions = []
    for record in read_bed_file(filename):
        if len(record) < 4:
            record.name = "%s*" % (record.contig,)

        regions.append(record)

    return regions


def parse_arguments(argv, ext):
    parser = argparse.ArgumentParser()

    parser.add_argument("infile", metavar="BAM",
                        help="Filename of a sorted BAM file. If set to '-' "
                             "the file is read from STDIN.")
    parser.add_argument("outfile", metavar="OUTPUT", nargs='?',
                        help="Filename of output table; defaults to name of "
                             "the input BAM with a '.depths' extension. If "
                             "set to '-' the table is printed to STDOUT.")
    parser.add_argument("--target-name", default=None,
                        help="Name used for 'Target' column; defaults to the "
                             "filename of the BAM file.")
    parser.add_argument("--regions-file", default=None, dest="regions_fpath",
                        help="BED file containing regions of interest; depth "
                             "is calculated only for these grouping by the "
                             "name used in the BED file, or the contig name "
                             "if no name has been specified for a record.")
    parser.add_argument('--max_contigs', default=100, type=int,
                        help="The maximum number of contigs allowed in a BAM "
                             "file. If this number is exceeded, the entire "
                             "set of contigs is aggregated into one pseudo-"
                             "contig named '<Genome>'. This is done to "
                             "limit table sizes [default: %(default)s]")
    parser.add_argument('--ignore-readgroups',
                        default=False, action="store_true",
                        help="Ignore readgroup information in reads, and only "
                             "provide aggreaged statistics; this is required "
                             "if readgroup information is missing or partial "
                             "[default: %(default)s]")

    args = parser.parse_args(argv)
    if not args.outfile:
        args.outfile = swap_ext(args.infile, ext)

    if args.ignore_readgroups:
        args.get_readgroup_func = _get_readgroup_ignored
    else:
        args.get_readgroup_func = _get_readgroup

    if not args.target_name:
        if args.infile == "-":
            args.target_name = "<STDIN>"
        else:
            args.target_name = os.path.basename(args.infile)

    return args


def main(process_func, argv, ext):
    # Set process name on supported platforms
    pypeline.common.system.set_procname()

    args = parse_arguments(argv, ext)
    args.regions = None
    if args.regions_fpath:
        args.regions = collect_bed_regions(args.regions_fpath)

    sys.stderr.write("Opening %r\n" % (args.infile,))
    with pysam.Samfile(args.infile) as handle:
        sort_bed_by_bamfile(handle, args.regions)
        return process_func(handle, args)


def _get_readgroup(record):
    try:
        return record.opt("RG")
    except KeyError:
        raise BAMStatsError("Record lacks readgroup tag ('RG'), use of "
                            "--ignore-readgroups is required to process this"
                            "file.\n    Record = %r" % (record.qname,))


def _get_readgroup_ignored(_):
    return None
