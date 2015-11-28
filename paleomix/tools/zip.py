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
Wrapper around gzip / bzip2.
This is motivated by the need to ensure proper handling of pipes, for which
the gzip / bzip2 commands do not behave consistently, to ensure that input
files are not removed, and finally to ensure that opening the pipe does not
block the main process.
"""
import os
import sys
import argparse
import subprocess


def main(argv):
    """Main function; takes a list of arguments but excluding sys.argv[0]."""
    prog = "paleomix zip"
    usage = "%s [options] in > out" % (prog,)

    parser = argparse.ArgumentParser(prog=prog, usage=usage)
    parser.add_argument("file", help="Input file.")
    parser.add_argument("--format", default="gz", choices=("bz2", "gz"),
                        help="Output format (gz / bz2).")

    args = parser.parse_args(argv)
    if args.format == "gz":
        command = ["gzip"]
    elif args.format == "bz2":
        command = ["bzip2"]
    else:
        assert False, args.format

    try:
        with open(args.file) as input_handle:
            subprocess.check_call(command,
                                  stdin=input_handle,
                                  preexec_fn=os.setsid,
                                  close_fds=True)

    except Exception, error:  # pylint: disable=W0703
        sys.stderr.write("Error running 'paleomix zip':\n    %s\n\n" % error)
        sys.stderr.write("Command = %s\n" % (" ".join(sys.argv),))
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
