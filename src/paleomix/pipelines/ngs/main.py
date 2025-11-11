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

import logging
import os
from argparse import Namespace

import paleomix.common.logging
import paleomix.node
from paleomix.common import resources
from paleomix.common.fileutils import swap_ext
from paleomix.common.yaml import YAMLError
from paleomix.pipeline import Pypeline
from paleomix.pipelines.ngs.config import build_parser
from paleomix.pipelines.ngs.pipeline import build_pipeline
from paleomix.pipelines.ngs.project import MakefileError, load_project


def main(argv: list[str]) -> int:
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


def _main_run(args: Namespace) -> int:
    if args.output is None:
        args.output = swap_ext(args.project, ".output")
    # FIXME: Add cli option
    args.temp_root = os.path.join(args.output, "cache", "temp")

    # FIXME: Place logs in output folder
    paleomix.common.logging.initialize_console_and_file_logging(
        log_level=args.log_level,
        log_color=args.log_color,
        log_file=args.log_file,
        auto_log_file=os.path.join(args.output, "cache", "logs", "pipeline"),
    )

    logger = logging.getLogger(__name__)

    try:
        logger.info("Reading project from %r", args.project)
        project = load_project(args.project)
    except (OSError, MakefileError, YAMLError) as error:
        logger.error("Error reading project: %s", error)
        return 1

    try:
        logger.info("Building pipeline for project")
        nodes = build_pipeline(args, project)
    except paleomix.node.NodeError as error:
        logger.error("Error while building pipeline: %s", error)
        return 1

    pipeline = Pypeline(
        nodes=nodes,
        group=args.work_group,
        temp_root=args.temp_root,
        max_threads=args.max_threads,
        intermediate_files=args.intermediate_files,
        required_files=args.require_files,
    )

    return pipeline.run(args.pipeline_mode)


def _main_new(_args: Namespace) -> int:
    print(resources.read_template("ngs.yaml"))

    return 0
