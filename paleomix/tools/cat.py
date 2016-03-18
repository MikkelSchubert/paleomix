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
Wrapper around cat / zcat / bzcat, which selects the appropriate commmand
based on the files specified on the command-line. Input files may be a mix
of different types of compressed / uncompressed files.
"""
import os
import sys
import argparse
import itertools
import subprocess


def _select_output(filename):
    """Returns a file-handle for 'filename'; if filename is '-' is stdout."""
    if filename in (None, '-'):
        return sys.stdout

    return open(filename, 'wb')


def _select_cat(filename):
    """Identifies the compression scheme of a given file (if any) and return a
    tuple with the appropriate cat command to decompress it.
    """
    with open(filename) as source:
        header = source.read(2)
        # The command "gzip -cd" is used instead of "zcat" because
        # OSX ships a broken zcat command (only accepts *.Z files).
        if header == "\x1f\x8b":
            return ("gzip", "-cd")
        elif header == "BZ":
            return ("bzip2", "-cd")
        return ("cat",)


def _call(input_files, output_file):
    """Call an appropriate cat on each input file, writing the contents to the
    file specified by 'output_file'; if the latter is '-', STDOUT is used.
    """
    with _select_output(output_file) as out_handle:
        for (command, filenames) in itertools.groupby(input_files,
                                                      _select_cat):
            command = list(command)
            command.extend(filenames)

            subprocess.check_call(command,
                                  stdout=out_handle,
                                  preexec_fn=os.setsid,
                                  close_fds=True)
    return 0


def parse_args(argv):
    parser = argparse.ArgumentParser(prog="paleomix cat")
    parser.add_argument("file", nargs="+",
                        help="One or more input files; these may be "
                             "uncompressed, compressed using gzip, or "
                             "compressed using bzip2.")
    parser.add_argument("--output", default=None,
                        help="Write output to this file; by default, output "
                             "is written to STDOUT.")

    return parser.parse_args(argv)


def main(argv):
    """Main function; takes a list of arguments but excluding sys.argv[0]."""
    args = parse_args(argv)

    try:
        return _call(input_files=args.file,
                     output_file=args.output)
    except Exception, error:
        sys.stderr.write("Error running 'paleomix cat':\n    %s\n\n" % error)
        sys.stderr.write("Command = %s\n" % (" ".join(sys.argv),))
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
