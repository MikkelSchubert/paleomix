#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import os
import sys
import time
import logging

import pypeline.ui
import pypeline.yaml
import pypeline.logger
import pypeline.tools.phylo_pipeline.makefile
import pypeline.tools.phylo_pipeline.mkfile as mkfile

from pypeline import Pypeline
from pypeline.common.console import print_err
from pypeline.tools.phylo_pipeline.makefile import \
     MakefileError, \
     read_makefiles
from pypeline.tools.phylo_pipeline.config import \
     ConfigError, \
     parse_config, \
     select_commands




def main(argv):
    try:
        config, args = parse_config(argv)
    except ConfigError, error:
        print_err(error)
        return 1

    if not args or ("help" in args):
        return 0
    elif (len(args) < 2) and ("mkfile" not in args):
        print_err("\nPlease specify at least one makefile!")
        return 1

    commands = select_commands(args.pop(0))
    if any((cmd == "mkfile") for (cmd, _) in commands):
        return mkfile.main(args[1:])

    if not os.path.exists(config.temp_root):
        try:
            os.makedirs(config.temp_root)
        except OSError, e:
            print_err("ERROR: Could not create temp root:\n\t%s" % (e,))
            return 1

    if not os.access(config.temp_root, os.R_OK | os.W_OK | os.X_OK):
        print_err("ERROR: Insufficient permissions for temp root: '%s'" % (config.temp_root,))
        return 1

    try:
        makefiles = read_makefiles(config, args, commands)
    except (MakefileError, pypeline.yaml.YAMLError, IOError), error:
        print_err("Error reading makefiles:",
                  "\n  %s:\n   " % (error.__class__.__name__,),
                  "\n    ".join(str(error).split("\n")),
                  file = sys.stderr)
        return 1

    logfile_template = time.strftime("phylo_pipeline.%Y%m%d_%H%M%S_%%02i.log")
    pypeline.logger.initialize(config, logfile_template)
    logger = logging.getLogger(__name__)

    pipeline = Pypeline(config)
    for (command_key, command_func) in commands:
        logger.info("Building %s pipeline ...", command_key)
        command_func(pipeline, config, makefiles)

    for makefile in makefiles:
        if "Nodes" in makefile:
            pipeline.add_nodes(makefile["Nodes"])

    if not pipeline.run(max_running = config.max_threads,
                        dry_run     = config.dry_run,
                        progress_ui = config.progress_ui):
        return 1
    return 0
