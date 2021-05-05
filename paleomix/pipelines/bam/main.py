import paleomix.resources

import paleomix.pipelines.bam.config as bam_config
import paleomix.pipelines.bam.pipeline as bam_pipeline


def main(argv, pipeline="bam"):
    if pipeline not in ("bam", "trim"):
        raise ValueError(pipeline)

    parser = bam_config.build_parser(pipeline)
    if not argv:
        parser.print_help()
        return 0

    args = parser.parse_args(argv)
    if args.command in ("new", "makefile", "mkfile"):
        return _main_template(args, pipeline=pipeline)
    elif args.command in ("example",):
        return paleomix.resources.copy_example("bam_pipeline", args)

    if args.command.startswith("dry"):
        args.dry_run = True

    return bam_pipeline.run(args, pipeline_variant=pipeline)


def _main_template(args, pipeline="bam"):
    if pipeline not in ("bam", "trim"):
        raise ValueError(pipeline)

    print(paleomix.resources.template("bam_head.yaml"))

    if pipeline == "bam":
        print(paleomix.resources.template("bam_options.yaml"))
        print()
        print(paleomix.resources.template("bam_prefixes.yaml"))

    print()
    print(paleomix.resources.template("bam_samples.yaml"))

    return 0
