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
from __future__ import print_function

import sys
import errno
import optparse
import fileinput

import pysam

import paleomix.common.vcffilter as vcffilter


def _read_files(filenames, args):
    in_header = True
    has_filters = False
    reset_filter = (args.reset_filter == 'yes')
    vcf_parser = pysam.asVCF()
    for line in fileinput.input(filenames):
        if not line.startswith("#"):
            in_header = False
            line = line.rstrip("\n\r")
            vcf = vcf_parser(line, len(line))
            if reset_filter:
                vcf.filter = '.'

            yield vcf
        elif in_header:
            if not (line.startswith("##") or has_filters):
                has_filters = True
                for item in sorted(vcffilter.describe_filters(args).items()):
                    print('##FILTER=<ID=%s,Description="%s">' % item)

            print(line, end="")


def main(argv):
    desc = "paleomix vcf_filter [options] [in1.vcf, ...]"
    parser = optparse.OptionParser(desc)
    parser.add_option('--reset-filter', default='no', choices=('yes', 'no'),
                      help="If set to 'yes', values in the 'FILTER' column "
                           "are cleared, and set according to the results "
                           "from running the filters implemented by this "
                           "tool. If set to 'no', any existing values are "
                           "retained, and any (new) failed filters are added "
                           "to these [default: %default].")

    vcffilter.add_varfilter_options(parser)
    (args, filenames) = parser.parse_args(argv)

    if (not filenames or "-" in filenames) and sys.stdin.isatty():
        parser.error("STDIN is a terminal, terminating!")

    try:
        for vcf in vcffilter.filter_vcfs(args, _read_files(filenames, args)):
            print(vcf)
    except IOError, error:
        # Check for broken pipe (head, less, etc).
        if error.errno != errno.EPIPE:
            raise


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
