#!/usr/bin/python3
#
# Copyright (c) 2013 Mikkel Schubert <MikkelSch@gmail.com>
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

import paleomix
import paleomix.common.logging

from paleomix.common.argparse import ArgumentParser, SUPPRESS
from paleomix.pipeline import add_io_argument_group, add_scheduling_argument_group

_DEFAULT_CONFIG_FILES = [
    "/etc/paleomix/phylo_pipeline.ini",
    "~/.paleomix/phylo_pipeline.ini",
]


def build_parser():
    parser = ArgumentParser(
        prog="paleomix phylo_pipeline",
        default_config_files=_DEFAULT_CONFIG_FILES,
    )
    parser.add_argument(
        "commands",
        help="One or more commands separated by '+'. Available commands are 'help', to "
        "display this message; 'example', to create an example project; 'new' to "
        "print a project template; 'genotype' to perform genotyping on a makefile; "
        "'msa' to perform multiple sequence alignment on a makefile; and 'phylogeny', "
        "to carry out phylogenetic inference on a makefile.",
    )
    parser.add_argument("files", nargs="*", help="One or more YAML files")

    paleomix.common.logging.add_argument_group(parser)

    group = add_scheduling_argument_group(parser)
    group.add_argument(
        "--examl-max-threads",
        default=1,
        type=int,
        help="Maximum number of threads for each instance of ExaML",
    )
    add_io_argument_group(parser)

    group = parser.add_argument_group("Required paths")
    group.add_argument(
        "--temp-root",
        default="./temp",
        type=os.path.abspath,
        help="Location for temporary files and folders",
    )
    group.add_argument(
        "--samples-root",
        default="./data/samples",
        help="Location of BAM files for each sample",
    )
    group.add_argument(
        "--regions-root",
        default="./data/regions",
        help="Location of BED files containing regions of interest",
    )
    group.add_argument(
        "--prefix-root",
        default="./data/prefixes",
        help="Location of prefixes (FASTAs)",
    )
    group.add_argument(
        "--destination",
        default="./results",
        help="The destination folder for result files",
    )

    # Removed options
    parser.add_argument("--refseq-root", help=SUPPRESS)
    parser.add_argument("--progress-ui", help=SUPPRESS)
    parser.add_argument("--ui-colors", help=SUPPRESS)

    return parser
