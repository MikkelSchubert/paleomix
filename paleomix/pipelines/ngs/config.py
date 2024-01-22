#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import multiprocessing
from enum import Enum

import paleomix
import paleomix.common.logging
import paleomix.pipeline
from paleomix.common.argparse import ArgumentParser, SubParsersAction

_DEFAULT_CONFIG_FILES = [
    "/etc/paleomix/ngs_pipeline.ini",
    "~/.paleomix/ngs_pipeline.ini",
]


class PipelineTarget(Enum):
    ALIGNMENTS = "alignments"
    HAPLOTYPES = "haplotypes"
    GENOTYPES = "genotypes"

    def __str__(self) -> str:
        return self.value


def build_parser() -> ArgumentParser:
    parser = ArgumentParser(prog="paleomix ngs")
    parser.set_defaults(command=None)

    subparsers = parser.add_subparsers()
    _build_run_parser(subparsers)
    _build_new_parser(subparsers)

    return parser


def _build_run_parser(subparsers: SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "run",
        help="Run pipeline on provided project file",
        default_config_files=_DEFAULT_CONFIG_FILES,
    )
    parser.set_defaults(command="run")

    parser.add_argument(
        "project",
        help="Run pipeline on these project files",
        metavar="FILE",
    )
    parser.add_argument(
        "output",
        nargs="?",
        metavar="FOLDER",
        help="Output folder for project files; defaults to the project file with an "
        "'.output' extension.",
    )

    parser.add_argument(
        "--target",
        metavar="OUTPUT",
        type=PipelineTarget,
        choices=list(PipelineTarget),
        default=PipelineTarget.GENOTYPES,
        help="Run the pipeline until the specified target output has been created",
    )

    parser.add_argument(
        "--bwa-algorithm",
        metavar="NAME",
        type=str.lower,
        choices=("mem", "mem2"),
        default="mem2",
        help="Mapping algorithm to use. Output should be identical, but BWA MEM2 is "
        "faster at the cost of 3x greater RAM usage",
    )

    paleomix.common.logging.add_argument_group(parser)

    group = paleomix.pipeline.add_scheduling_argument_group(parser)

    threaded_software = {
        "fastp": 4,
        "BWA": 8,
        "SAMtools": 3,
    }

    for name, max_default in threaded_software.items():
        group.add_argument(
            f"--max-threads-{name.lower()}",
            type=int,
            default=max(1, min(max_default, multiprocessing.cpu_count() // 2)),
            help=f"Max number of threads to use per {name} instance",
        )

    paleomix.pipeline.add_io_argument_group(parser)


def _build_new_parser(subparsers: SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser("new")
    parser.set_defaults(command="new")
