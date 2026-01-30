# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import multiprocessing
import os

import paleomix
import paleomix.common.logging
import paleomix.pipeline
from paleomix.common.argparse import SUPPRESS, ArgumentParser, SubParsersAction

_DEFAULT_CONFIG_FILES = [
    "/etc/paleomix/bam_pipeline.ini",
    "~/.paleomix/bam_pipeline.ini",
]


def build_parser(pipeline_variant: str) -> ArgumentParser:
    parser = ArgumentParser(prog=f"paleomix {pipeline_variant}")

    subparsers = parser.add_subparsers(dest="command", metavar="command")
    add_makefile_command(subparsers)
    add_run_command(subparsers)
    add_copy_example_command(subparsers)

    return parser


def add_makefile_command(subparsers: SubParsersAction[ArgumentParser]) -> None:
    subparsers.add_parser(
        "new",
        help="Print project template",
        aliases=("mkfile", "makefile"),
    )


def add_run_command(subparsers: SubParsersAction[ArgumentParser]) -> None:
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

    # Removed options
    group.add_argument("--jar-root", help=SUPPRESS)
    parser.add_argument("--jre-option", "--jre-options", help=SUPPRESS)
    parser.add_argument("--gatk-max-threads", help=SUPPRESS)
    parser.add_argument("--progress-ui", help=SUPPRESS)
    parser.add_argument("--ui-colors", help=SUPPRESS)


def add_copy_example_command(subparsers: SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser("example", help="Create example project")

    parser.add_argument(
        "destination",
        default=".",
        nargs="?",
        help="Destination folder for example data.",
    )
