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
import optparse

import yaml

import pypeline.ui
import pypeline.logger
import pypeline.tools.phylo_pipeline.makefile
import pypeline.tools.phylo_pipeline.mkfile as mkfile
import pypeline.tools.phylo_pipeline.parts.genotype as genotype
import pypeline.tools.phylo_pipeline.parts.msa as msa
import pypeline.tools.phylo_pipeline.parts.paml as paml
import pypeline.tools.phylo_pipeline.parts.phylo as phylo

from pypeline import Pypeline
from pypeline.common.console import print_err
from pypeline.tools.phylo_pipeline.makefile import \
     read_makefiles


_COMMANDS = {
    "mkfile"          : True,
    "genotype"        : genotype.chain,
    "genotyping"      : genotype.chain,
    "msa"             : msa.chain,
    "paml:codeml"     : paml.chain_codeml,
    "phylogeny:examl" : phylo.chain_examl,
    }


class CustomHelpFormatter(optparse.IndentedHelpFormatter):
    def format_description(self, description):
        return description or ""


def build_options_parser():
    parser = optparse.OptionParser("%prog <command> [options] [makefiles]")
    parser.formatter = CustomHelpFormatter()
    parser.formatter.set_parser(parser)
    parser.description = \
      "Commands:\n" \
      "  -- %prog help            -- Display this message.\n" \
      "  -- %prog mkfile [...]    -- Print makefile template.\n" \
      "  -- %prog genotype [...]  -- Carry out genotyping according to makefile.\n" \
      "  -- %prog msa [...]       -- Carry out multiple sequence alignment.\n" \
      "  -- %prog phylogeny [...] -- Carry out phylogenetic inference.\n"

    pypeline.ui.add_optiongroup(parser)
    pypeline.logger.add_optiongroup(parser)

    group  = optparse.OptionGroup(parser, "Scheduling")
    group.add_option("--max-threads",        default = 12, type = int,
                     help = "Maximum number of threads to use in total [%default]")
    group.add_option("--max-examl-threads",  default = 4, type = int,
                     help = "Maximum number of threads to use for each instance of ExaML [%default]")
    group.add_option("--dry-run",            default = False, action="store_true",
                     help = "If passed, only a dry-run in performed, the dependency tree is printed, and no tasks are executed.")
    parser.add_option_group(group)

    group  = optparse.OptionGroup(parser, "Required paths")
    group.add_option("--temp-root",          default = "./temp",
                     help = "Location for temporary files and folders [%default]")
    group.add_option("--samples-root",       default = "./data/samples",
                     help = "Location of BAM files for each sample.")
    group.add_option("--intervals-root",     default = "./data/intervals",
                     help = "Location of BED files containing intervals of interest [%default]")
    group.add_option("--genomes-root",       default = "./data/genomes",
                     help = "Location of reference genomes (FASTAs) [%default]")
    group.add_option("--destination",        default = "./results",
                     help = "The destination folder for result files [%default]")
    parser.add_option_group(group)

    return parser


def parse_args(argv):
    options_parser = build_options_parser()
    options, args  = options_parser.parse_args(argv)

    if (len(args) < 2) and (args != ["mkfile"]):
        print_err("Please specify at least one analysis step and one makefile!", file = sys.stderr)
        return None, None

    commands = _select_commands(args[0])
    if any((func is None) for (_, func) in commands):
        unknown_commands = ", ".join(repr(key) for (key, func) in commands if func is None)
        print_err("Unknown analysis step(s): %s" % (unknown_commands,), file = sys.stderr)
        return None, None

    if not os.path.exists(options.temp_root):
        try:
            os.makedirs(options.temp_root)
        except OSError, e:
            print_err("ERROR: Could not create temp root:\n\t%s" % (e,), file = sys.stderr)
            return None, None

    if not os.access(options.temp_root, os.R_OK | os.W_OK | os.X_OK):
        print_err("ERROR: Insufficient permissions for temp root: '%s'" \
                  % (options.temp_root,), file = sys.stderr)
        return None, None

    return options, args



def _select_commands(chain):
    commands = []
    for command in chain.split("+"):
        command_key  = command.strip().lower()
        command_func = None

        if command in _COMMANDS:
            command_func = _COMMANDS[command_key]
        elif len(command) >= 3:
            for (key, value) in _COMMANDS.iteritems():
                if key.startswith(command):
                    command_key  = key
                    command_func = value
                    break

        commands.append((command_key, command_func))

    return commands


def main(argv):
    options, args = parse_args(argv)
    if not (args and options) or (args and (args[0] == "help")):
        return 1

    commands = _select_commands(args.pop(0))
    if any((cmd == "mkfile") for (cmd, _) in commands):
        return mkfile.main(args[1:])

    try:
        makefiles = read_makefiles(args)
    except (StandardError, yaml.YAMLError), error:
        print_err("Error reading makefiles:\n    ",
                  "\n    ".join(str(error).split("\n")),
                  file = sys.stderr)
        return 1

    logfile_template = time.strftime("phylo_pipeline.%Y%m%d_%H%M%S_%%02i.log")
    pypeline.logger.initialize(options, logfile_template)
    logger = logging.getLogger(__name__)

    pipeline = Pypeline(options)
    for (command_key, command_func) in commands:
        logger.info("Building %s pipeline ...", command_key)
        command_func(pipeline, options, makefiles)

    for makefile in makefiles:
        if "Nodes" in makefile:
            pipeline.add_nodes(makefile["Nodes"])

    if not pipeline.run(max_running = options.max_threads,
                        dry_run     = options.dry_run,
                        progress_ui = options.progress_ui):
        return 1
    return 0
