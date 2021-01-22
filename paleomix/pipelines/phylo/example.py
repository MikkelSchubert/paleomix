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
import os
import sys

import paleomix.common.logging
import paleomix.resources
import paleomix.yaml


def build_parser():
    parser = argparse.ArgumentParser(prog="paleomix phylo_pipeline")

    parser.add_argument(
        "destination",
        default=".",
        nargs="?",
        help="Destination folder for example data.",
    )

    paleomix.common.logging.add_argument_group(parser)

    return parser


def main(argv):
    parser = build_parser()
    config = parser.parse_args(argv)
    paleomix.common.logging.initialize_console_logging()

    if paleomix.resources.copy_example("phylo_pipeline", config):
        return 1

    # Update interpreter to match the one currently in use;
    # this is required since we may be running from a virtual env
    filename = os.path.join(config.destination, "phylo_pipeline", "synthesize_reads.py")

    with open(filename) as handle:
        header, lines = handle.read().split("\n", 1)

    with open(filename, "w") as handle:
        handle.write("#!%s\n" % (os.path.abspath(sys.executable)))
        handle.write(lines)

    return 0
