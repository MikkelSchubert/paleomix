import multiprocessing

import paleomix
import paleomix.common.logging

from paleomix.common.argparse import ArgumentParser


_DEFAULT_CONFIG_FILES = [
    "/etc/paleomix/ngs_pipeline.ini",
    "~/.paleomix/ngs_pipeline.ini",
]


def build_parser():
    parser = ArgumentParser(prog="paleomix ngs")
    parser.set_defaults(command=None)

    subparsers = parser.add_subparsers()
    _build_run_parser(subparsers)
    _build_new_parser(subparsers)

    return parser


def _build_run_parser(subparsers):
    parser = subparsers.add_parser("run")
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

    paleomix.common.logging.add_argument_group(parser)

    group = parser.add_argument_group("Scheduling")
    group.add_argument(
        "--dry-run",
        action="store_true",
        help="Build pipeline and check prerequisites, but do not execute any tasks",
    )
    group.add_argument(
        "--max-threads",
        type=int,
        default=max(2, multiprocessing.cpu_count()),
        help="Max number of threads to use in total",
    )

    group = parser.add_argument_group("Input/Output")
    group.add_argument(
        "--list-output-files",
        action="store_true",
        default=False,
        help="List all output files generated by pipeline for the makefile(s).",
    )


def _build_new_parser(subparsers):
    parser = subparsers.add_parser("new")
    parser.set_defaults(command="new")
