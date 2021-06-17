#!/usr/bin/env python3
# -*- coding: utf8 -*-
import sys
from pathlib import Path

from paleomix.common.argparse import ArgumentParser
from paleomix.common.formats.bed import (
    merge_bed_records,
    pad_bed_records,
    read_bed_file,
)


def main_pad(args):
    contigs = {}
    with args.fai.open("rt") as handle:
        for line in handle:
            name, length, _ = line.split("\t", 2)
            if name in contigs:
                print(
                    "Reference genome contains multiple identically named contigs: %r"
                    % (name,),
                    file=sys.stderr,
                )

            contigs[name] = int(length)

    records = []
    for filepath in args.bed:
        records.extend(read_bed_file(filepath, contigs=contigs))

    padded_records = pad_bed_records(
        records=records,
        padding=args.padding,
        max_sizes=contigs,
    )

    for record in merge_bed_records(padded_records):
        print(record)

    return 0


def build_parser():
    parser = ArgumentParser("paleomix :bedtools")
    parser.set_defaults(main=None)

    subparsers = parser.add_subparsers()

    subparser = subparsers.add_parser("pad")
    subparser.set_defaults(main=main_pad)
    subparser.add_argument("fai", type=Path)
    subparser.add_argument("bed", nargs="+", type=Path)
    subparser.add_argument("--padding", default=0, type=int)

    return parser


def main(argv):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.main is None:
        parser.print_usage()
        return 0

    return args.main(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
