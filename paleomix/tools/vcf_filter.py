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
import errno
import sys

import pysam

import paleomix.common.argparse as argparse
import paleomix.common.vcffilter as vcffilter

from paleomix.common.fileutils import open_ro


def _read_files(args):
    in_header = True
    has_filters = False
    vcf_parser = pysam.asVCF()
    for filename in args.filenames:
        with open_ro(filename, "rb") as handle:
            for line in handle:
                if not line.startswith(b"#"):
                    in_header = False
                    line = line.rstrip(b"\n\r")
                    vcf = vcf_parser(line, len(line))
                    if args.reset_filter:
                        vcf.filter = "."

                    yield vcf
                elif in_header:
                    if not (line.startswith(b"##") or has_filters):
                        has_filters = True
                        for item in sorted(vcffilter.describe_filters(args).items()):
                            print('##FILTER=<ID=%s,Description="%s">' % item)

                    print(line.decode("utf-8"), end="")


def main(argv):
    parser = argparse.ArgumentParser(prog="paleomix vcf_filter")

    parser.add_argument(
        "filenames",
        nargs="*",
        help="VCF files; may be gzip/bzip2 compresssed. Leave blank or use '-' to "
        "read from STDIN",
        metavar="file",
    )

    parser.add_argument(
        "--reset-filter",
        default=False,
        action="store_true",
        help="If set, values in the 'FILTER' column  are cleared, and set according "
        "to the results from running the filters implemented by this tool. If not "
        "set, any existing values are retained, and any (new) failed filters are "
        "added to these.",
    )

    vcffilter.add_varfilter_options(parser)
    args = parser.parse_args(argv)

    if (not args.filenames or "-" in args.filenames) and sys.stdin.isatty():
        parser.error("STDIN is a terminal, terminating!")

    try:
        for vcf in vcffilter.filter_vcfs(args, _read_files(args)):
            print(vcf)
    except IOError as error:
        # Check for broken pipe (head, less, etc).
        if error.errno != errno.EPIPE:
            raise


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
