#!/usr/bin/python3
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

import paleomix.common.logging
import paleomix.pipelines.phylo.example as example
import paleomix.pipelines.phylo.mkfile as mkfile
import paleomix.pipelines.phylo.parts.genotype as genotype
import paleomix.pipelines.phylo.parts.msa as msa
import paleomix.pipelines.phylo.parts.phylo as phylo
import paleomix.yaml

from paleomix.pipeline import Pypeline
from paleomix.pipelines.phylo.config import build_parser
from paleomix.pipelines.phylo.makefile import MakefileError, read_makefiles


_COMMANDS = {
    "genotype": genotype.chain,
    "msa": msa.chain,
    "phylogeny": phylo.chain_examl,
}


def main(argv):
    log = logging.getLogger(__name__)

    commands = [it.lower().strip() for it in argv[0].split("+")] if argv else []
    if "example" in commands:
        return example.main(argv[1:])

    parser = build_parser()
    config = parser.parse_args(argv)
    config.commands = commands
    if "help" in config.commands:
        parser.print_help()
        return 0
    elif "makefile" in config.commands or "mkfile" in config.commands:
        return mkfile.main(config.files)

    commands = []
    for key in config.commands:
        func = _COMMANDS.get(key)
        if func is None:
            log.error("unknown command %r", key)
            return 1

        commands.append((key, func))

    paleomix.common.logging.initialize(
        log_level=config.log_level,
        log_file=config.log_file,
        name="phylo_pipeline",
    )

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
        makefiles = read_makefiles(config, commands)
    except (MakefileError, paleomix.yaml.YAMLError, IOError) as error:
        log.error("Error reading makefiles:\n%s", error)
        return 1

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
