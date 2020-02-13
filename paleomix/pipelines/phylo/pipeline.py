#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import logging
import os
import sys

import paleomix.logger
import paleomix.resources
import paleomix.pipelines.phylo.mkfile as mkfile
import paleomix.yaml
from paleomix.pipeline import Pypeline
from paleomix.pipelines.phylo.config import (
    ConfigError,
    parse_config,
    select_commands,
)
from paleomix.pipelines.phylo.makefile import MakefileError, read_makefiles


def main(argv):
    log = logging.getLogger(__name__)

    try:
        config, args = parse_config(argv)
    except ConfigError as error:
        log.error("%s", error)
        return 1

    if not args or ("help" in args):
        return 0
    elif args[0] in ("example", "examples"):
        if paleomix.resources.copy_example("phylo_pipeline", argv[1:]):
            return 1

        # Update interpreter to match the one currently in use;
        # this is required since we may be running from a virtual env
        filename = os.path.join(argv[1], "phylo_pipeline", "synthesize_reads.py")

        with open(filename) as handle:
            header, lines = handle.read().split("\n", 1)

        with open(filename, "w") as handle:
            handle.write("#!%s\n" % (os.path.abspath(sys.executable)))
            handle.write(lines)

        return 0
    elif (len(args) < 2) and ("mkfile" not in args and "makefile" not in args):
        log.error("\nPlease specify at least one makefile!")
        return 1

    commands = select_commands(args.pop(0))
    if any((cmd in ("makefile", "mkfile")) for (cmd, _) in commands):
        return mkfile.main(args[1:])

    if not os.path.exists(config.temp_root):
        try:
            os.makedirs(config.temp_root)
        except OSError as error:
            log.error("Could not create temp root:\n\t%s", error)
            return 1

    if not os.access(config.temp_root, os.R_OK | os.W_OK | os.X_OK):
        log.error("Insufficient permissions for temp root: %r", config.temp_root)
        return 1

    # Init worker-threads before reading in any more data
    pipeline = Pypeline(config)

    try:
        makefiles = read_makefiles(config, args, commands)
    except (MakefileError, paleomix.yaml.YAMLError, IOError) as error:
        log.error("Error reading makefiles:\n%s", error)
        return 1

    paleomix.logger.initialize(
        log_level=config.log_level, log_file=config.log_file, name="phylo_pipeline"
    )

    for (command_key, command_func) in commands:
        log.info("Building %s pipeline", command_key)
        command_func(pipeline, config, makefiles)

    for makefile in makefiles:
        if "Nodes" in makefile:
            pipeline.add_nodes(makefile["Nodes"])

    if config.list_input_files:
        log.info("Printing output files")
        pipeline.print_input_files()
        return 0
    elif config.list_output_files:
        log.info("Printing output files")
        pipeline.print_output_files()
        return 0
    elif config.list_executables:
        log.info("Printing required executables")
        pipeline.print_required_executables()
        return 0

    if not pipeline.run(max_threads=config.max_threads, dry_run=config.dry_run):
        return 1
    return 0
