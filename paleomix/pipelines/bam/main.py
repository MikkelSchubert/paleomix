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

from typing import List

import paleomix.pipelines.bam.config as bam_config
import paleomix.pipelines.bam.pipeline as bam_pipeline
import paleomix.resources


def main(argv: List[str], pipeline: str = "bam") -> int:
    if pipeline not in ("bam", "trim"):
        raise ValueError(pipeline)

    parser = bam_config.build_parser(pipeline)
    if not argv:
        parser.print_help()
        return 0

    args = parser.parse_args(argv)
    if args.command in ("new", "makefile", "mkfile"):
        return _main_template()
    elif args.command in ("example",):
        return paleomix.resources.copy_example("bam_pipeline", args.destination)

    if args.command.startswith("dry") and args.pipeline_mode == "run":
        args.pipeline_mode = "dry_run"

    args.pipeline_variant = pipeline

    return bam_pipeline.run(args)


def _main_template() -> int:
    print(paleomix.resources.template("bam_head.yaml"))
    print(paleomix.resources.template("bam_options.yaml"))

    print()

    print(paleomix.resources.template("bam_prefixes.yaml"))

    print()
    print(paleomix.resources.template("bam_samples.yaml"))

    return 0
