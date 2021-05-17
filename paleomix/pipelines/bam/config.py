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
import multiprocessing

import paleomix
import paleomix.pipeline
import paleomix.common.logging

from paleomix.resources import add_copy_example_command
from paleomix.common.argparse import ArgumentParser, SUPPRESS


_DEFAULT_CONFIG_FILES = [
    "/etc/paleomix/bam_pipeline.ini",
    "~/.paleomix/bam_pipeline.ini",
]


def build_parser(pipeline_variant):
    parser = ArgumentParser(prog="paleomix %s" % (pipeline_variant,))

    subparsers = parser.add_subparsers(dest="command", metavar="command")
    add_makefile_command(subparsers)
    add_run_command(subparsers)
    add_copy_example_command(subparsers)

    return parser


def add_makefile_command(subparsers):
    subparsers.add_parser(
        "new",
        help="Print project template",
        aliases=("mkfile", "makefile"),
    )


def add_run_command(subparsers):
    parser = subparsers.add_parser(
        "run",
        aliases=("dryrun",),
        help="Run pipeline on provided makefiles",
        default_config_files=_DEFAULT_CONFIG_FILES,
    )

    parser.add_argument(
        "makefiles",
        nargs="+",
        help="Run pipeline on these makefiles",
        metavar="makefile",
    )

    paleomix.common.logging.add_argument_group(parser)

    group = paleomix.pipeline.add_scheduling_argument_group(parser)
    group.add_argument(
        "--adapterremoval-max-threads",
        type=int,
        default=min(3, multiprocessing.cpu_count()),
        help="Max number of threads to use per AdapterRemoval instance",
    )
    group.add_argument(
        "--bowtie2-max-threads",
        type=int,
        default=max(1, min(8, multiprocessing.cpu_count() // 2)),
        help="Max number of threads to use per Bowtie2 instance",
    )
    group.add_argument(
        "--bwa-max-threads",
        type=int,
        default=max(1, min(8, multiprocessing.cpu_count() // 2)),
        help="Max number of threads to use per BWA instance",
    )
    group.add_argument(
        "--samtools-max-threads",
        type=int,
        default=max(1, min(3, multiprocessing.cpu_count() // 2)),
        help="Max number of threads to use per SAMTools instance",
    )

    paleomix.pipeline.add_io_argument_group(parser)

    group = parser.add_argument_group("Required paths")
    group.add_argument(
        "--jar-root",
        default=os.path.expanduser("~/install/jar_root"),
        help="Folder containing Picard JARs (http://picard.sf.net)",
    )
    group.add_argument(
        "--temp-root",
        default="./temp/",
        type=os.path.abspath,
        help="Location for temporary files and folders",
    )
    group.add_argument(
        "--destination",
        default=".",
        help="The destination folder for result files.",
    )

    group = parser.add_argument_group("Misc")
    group.add_argument(
        "--jre-option",
        "--jre-options",
        metavar="OPTION",
        dest="jre_options",
        action="append",
        default=[],
        help="May be specified one or more times with options to be passed "
        "to the JRE (Jave Runtime Environment); e.g. to change the "
        "maximum amount of memory (default is -Xmx4g)",
    )

    # Removed options
    parser.add_argument("--gatk-max-threads", help=SUPPRESS)
    parser.add_argument("--progress-ui", help=SUPPRESS)
    parser.add_argument("--ui-colors", help=SUPPRESS)
