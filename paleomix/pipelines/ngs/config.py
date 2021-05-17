import paleomix
import paleomix.pipeline
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
    paleomix.pipeline.add_argument_groups(parser)


def _build_new_parser(subparsers):
    parser = subparsers.add_parser("new")
    parser.set_defaults(command="new")
