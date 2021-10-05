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

    return bam_pipeline.run(args, pipeline_variant=pipeline)


def _main_template() -> int:
    print(paleomix.resources.template("bam_head.yaml"))
    print(paleomix.resources.template("bam_options.yaml"))

    print()

    print(paleomix.resources.template("bam_prefixes.yaml"))

    print()
    print(paleomix.resources.template("bam_samples.yaml"))

    return 0
