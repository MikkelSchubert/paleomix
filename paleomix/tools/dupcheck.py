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


import argparse
import sys

import paleomix.nodes.validation as validation


class ErrHandler:
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


def parse_args(argv):
    parser = argparse.ArgumentParser(prog="paleomix dupcheck")
    parser.add_argument("files", nargs="+", help="One or more input BAM files.")
    parser.add_argument(
        "--quiet",
        default=False,
        action="store_true",
        help="Only print the number of BAM records where 1 or "
        "more potential duplicates duplicates were "
        "identified.",
    )

    return parser.parse_args(argv)


def main(argv):
    """Main function; takes a list of arguments but excluding sys.argv[0]."""
    args = parse_args(argv)
    handler = ErrHandler(quiet=args.quiet)
    validation.check_bam_files(args.files, handler)

    if args.quiet:
        print("%i" % (handler.duplicate_reads,))
    else:
        print("Found %i record(s) with duplicates." % (handler.duplicate_reads,))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
