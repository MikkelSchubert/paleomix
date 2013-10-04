#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import optparse

import pypeline

import pypeline.tools.phylo_pipeline.parts.genotype as genotype
import pypeline.tools.phylo_pipeline.parts.msa as msa
import pypeline.tools.phylo_pipeline.parts.paml as paml
import pypeline.tools.phylo_pipeline.parts.phylo as phylo

from pypeline.config import \
     ConfigError, \
     PerHostValue, \
     PerHostConfig




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


def select_commands(chain):
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


def _run_config_parser(argv):
    per_host_cfg = PerHostConfig("phylo_pipeline")
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
    group.add_option("--examl-max-threads",  default = PerHostValue(4), type = int,
                     help = "Maximum number of threads to use for each instance of ExaML [%default]")
    group.add_option("--max-threads",        default = per_host_cfg.max_threads, type = int,
                     help = "Maximum number of threads to use in total [%default]")
    group.add_option("--dry-run",            default = False, action="store_true",
                     help = "If passed, only a dry-run in performed, the dependency tree is printed, "
                            "and no tasks are executed.")
    parser.add_option_group(group)

    group  = optparse.OptionGroup(parser, "Required paths")
    group.add_option("--temp-root",    default = per_host_cfg.temp_root,
                     help = "Location for temporary files and folders [%default]")
    group.add_option("--samples-root", default = PerHostValue("./data/samples", is_path = True),
                     help = "Location of BAM files for each sample [%default]")
    group.add_option("--regions-root", default = PerHostValue("./data/regions", is_path = True),
                     help = "Location of BED files containing regions of interest [%default]")
    group.add_option("--prefix-root",  default = PerHostValue("./data/prefixes", is_path = True),
                     help = "Location of prefixes (FASTAs) [%default]")
    group.add_option("--destination",  default = "./results",
                     help = "The destination folder for result files [%default]")
    parser.add_option_group(group)

    return per_host_cfg.parse_args(parser, argv)



def parse_config(argv):
    options, args  = _run_config_parser(argv)
    print options.max_threads

    if (len(args) < 2) and (args != ["mkfile"]):
        raise ConfigError("Please specify at least one analysis step and one makefile!")

    commands = select_commands(args[0])
    if any((func is None) for (_, func) in commands):
        unknown_commands = ", ".join(repr(key) for (key, func) in commands if func is None)
        raise ConfigError("Unknown analysis step(s): %s" % (unknown_commands,))

    if not os.path.exists(options.temp_root):
        try:
            os.makedirs(options.temp_root)
        except OSError, e:
            raise ConfigError("ERROR: Could not create temp root:\n\t%s" % (e,))

    if not os.access(options.temp_root, os.R_OK | os.W_OK | os.X_OK):
        raise ConfigError("ERROR: Insufficient permissions for temp root: '%s'" \
                          % (options.temp_root,))

    return options, args
