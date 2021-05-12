import os

import logging

import paleomix.common.logging

from paleomix.common.fileutils import swap_ext

from paleomix.pipeline import Pypeline
from paleomix.pipelines.ngs.config import build_parser
from paleomix.pipelines.ngs.project import load_project, MakefileError
from paleomix.pipelines.ngs.pipeline import build_pipeline
from paleomix import resources


def main(argv):
    parser = build_parser()
    if not argv:
        parser.print_help()
        return 0

    args = parser.parse_args(argv)
    if args.command == "run":
        return _main_run(args)
    elif args.command == "new":
        return _main_new(args)

    parser.print_help()
    return 0


def _main_run(args):
    if args.output is None:
        args.output = swap_ext(args.project, ".output")
    # FIXME: Add cli option
    args.temp_root = os.path.join(args.output, "cache", "temp")

    # FIXME: Place logs in output folder
    paleomix.common.logging.initialize(
        log_level=args.log_level,
        log_file=args.log_file,
        auto_log_file=os.path.join(args.output, "cache", "logs", "pipeline"),
    )

    logger = logging.getLogger(__name__)

    # Initialize worker-threads before reading in any more data
    pipeline = Pypeline(args)

    try:
        logger.info("Reading project from %r", args.project)
        project = load_project(args.project)
    except (MakefileError, paleomix.yaml.YAMLError, IOError) as error:
        logger.error("Error reading project: %s", error)
        return 1

    if not project["Samples"]:
        logger.warning(
            "Project does not contain any samples; genomes will be prepared for later "
            "use, but no reads mapped to these"
        )

    try:
        logger.info("Building pipeline for project")
        nodes = build_pipeline(args, project)
    except paleomix.node.NodeError as error:
        logger.error("Error while building pipeline: %s", error)
        return 1

    pipeline.add_nodes(nodes)

    if args.list_output_files:
        logger.info("Printing output files")
        pipeline.print_output_files()
        return 0

    logger.info("Running pipeline")
    return not pipeline.run(dry_run=args.dry_run, max_threads=args.max_threads)


def _main_new(args):
    print(resources.template("ngs.yaml"))

    return 0
