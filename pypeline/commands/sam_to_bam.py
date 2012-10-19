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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
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
An alternative to 'samtools view -Sb -', needed because samtools does not
properly handle SAM-files with no records. For example, the following
command will not work:
$ samtools view -H INPUT | samtools view -Sb -
"""

import os
import sys
import optparse

import pysam

import pypeline.ui
import pypeline.atomiccmd
from pypeline.common.utilities import set_in


def call(min_quality, exclude_flags, flag_as_sorted):
    input_file = pysam.Samfile("-", "r")
    header = input_file.header
    if flag_as_sorted:
        set_in(header, ("HD", "SO"), "coordinate")
        set_in(header, ("HD", "GO"), "none")
        set_in(header, ("HD", "VN"), "1.0")

    output_file = pysam.Samfile("-", "wb", header = header)
    for read in input_file:
        if read.flag & exclude_flags:
            continue
        elif read.mapq < min_quality:
            continue
               
        output_file.write(read)
    output_file.close()
    input_file.close()


def main(argv):
    parser = optparse.OptionParser()
    parser.add_option("-q", "--min-quality", type = int, default = 0)
    parser.add_option("-F", "--exclude-flags", type = int, default = 0)
    parser.add_option("--flag-as-sorted", default = False, action = "store_true")

    config, args = parser.parse_args(argv)
    if args:
        pypeline.ui.print_err("%s does not take any arguments." % sys.argv[0])
        sys.exit(1)
                               
    return call(min_quality    = config.min_quality,
                exclude_flags  = config.exclude_flags,
                flag_as_sorted = config.flag_as_sorted)


def build_atomiccmd(cls, min_quality, exclude_flags, stdin, output_file, flag_as_sorted = False):
    filename = os.path.realpath(__file__)
    if filename.endswith(".pyc"):
        filename = filename[:-1]

    call = [filename, 
            "-q", min_quality,
            "-F", exclude_flags]
    if flag_as_sorted:
        call.append("--flag-as-sorted")

    return cls(call,
               IN_STDIN   = stdin,
               OUT_STDOUT = output_file)
     

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
 