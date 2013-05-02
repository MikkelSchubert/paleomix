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
import pwd
import sys
import copy
import string
import optparse
import collections
import ConfigParser

import pypeline
import pypeline.ui as ui

from pypeline.node import MetaNode
from pypeline.nodes.picard import BuildSequenceDictNode
from pypeline.nodes.samtools import FastaIndexNode, BAMIndexNode

from pypeline.common.fileutils import swap_ext, add_postfix

from pypeline.tools.bam_pipeline.makefile import *
from pypeline.tools.bam_pipeline.nodes import MapDamageNode

import pypeline.tools.bam_pipeline.parts as parts




def _add_mapdamage_nodes(config, makefile, target):
    if "mapDamage" not in makefile["Options"]["Features"]:
        return

    nodes = []
    for prefix in target.prefixes:
        libraries = []
        for sample in prefix.samples:
            for library in sample.libraries:
                folder = os.path.join(config.destination, "%s.%s.mapDamage" % (target.name, prefix.name), library.name)
                libraries.append(MapDamageNode(config           = config,
                                               reference        = prefix.reference,
                                               input_files      = library.bams.keys(),
                                               output_directory = folder,
                                               dependencies     = library.bams.values()))

        nodes.append(MetaNode(description = prefix.name,
                              subnodes    = libraries))

    target.add_extra_nodes("mapDamage", nodes)



def _add_extra_nodes(config, makefile, targets):
    for target in targets:
        _add_mapdamage_nodes(config, makefile, target)
        parts.add_statistics_nodes(config, makefile, target)

    return targets


def build_pipeline_trimming(config, makefile):
    """Builds only the nodes required to produce trimmed reads.
    This reduces the required complexity of the makefile to a minimum."""

    nodes = []
    for prefix in makefile["Prefixes"].itervalues():
        for (target, samples) in makefile["Targets"].iteritems():
            for (sample, libraries) in samples.iteritems():
                for (library, barcodes) in libraries.iteritems():
                    for (barcode, record) in barcodes.iteritems():
                        lane = parts.Lane(config, prefix, record, barcode)
                        if lane.reads and lane.reads.nodes:
                            nodes.extend(lane.reads.nodes)
        break # Only one prefix is required
    return nodes


def build_pipeline_full(config, makefile, return_nodes = True):
    targets = []
    features = makefile["Options"]["Features"]
    for (target_name, sample_records) in makefile["Targets"].iteritems():
        prefixes = []
        for (prefix_name, prefix) in makefile["Prefixes"].iteritems():
            samples = []
            for (sample_name, library_records) in sample_records.iteritems():
                libraries = []
                for (library_name, barcode_records) in library_records.iteritems():
                    lanes = []
                    for (barcode, record) in barcode_records.iteritems():
                        lanes.append(parts.Lane(config, prefix, record, barcode))

                    if any(lane.bams for lane in lanes):
                        libraries.append(parts.Library(config, prefix, lanes, library_name))

                if libraries:
                    samples.append(parts.Sample(config, prefix, libraries, sample_name))

            if samples:
                prefixes.append(parts.Prefix(config, prefix, samples, features, target_name))

        if prefixes:
            targets.append(parts.Target(config, prefixes, target_name))

    targets = _add_extra_nodes(config, makefile, targets)
    if not return_nodes:
        return targets

    return [target.node for target in targets]



def _make_target_list(config, makefiles):
    target_list = {}
    for target in build_pipeline_full(config, makefiles, return_nodes = False):
        target_list[(target.name,)] = [target.node]

        for prefix in target.prefixes:
            target_list[(target.name, prefix.name)] = prefix.bams.values()

            for sample in prefix.samples:
                target_list[(target.name, prefix.name, sample.name)] = sample.bams.values()

                for library in sample.libraries:
                    target_list[(target.name, prefix.name, sample.name, library.name)] = library.bams.values()

                    for lane in library.lanes:
                        lane_bams = []
                        for files_and_nodes in lane.bams.itervalues():
                            lane_bams.extend(files_and_nodes.itervalues())

                        target_list[(target.name, prefix.name, sample.name, library.name, lane.name)] = lane_bams

                        for (reads_type, bams) in lane.bams.iteritems():
                            target_list[(target.name, prefix.name, sample.name, library.name, lane.name, reads_type)] = bams.values()

                        if lane.reads and lane.reads.nodes:
                            target_list[(target.name, "reads", sample.name, library.name, lane.name)] = lane.reads.nodes

    return target_list


def list_targets_for(config, makefiles, show):
    target_list = _make_target_list(config, makefiles)
    length = {"targets"   : 1, "prefixes" : 2, "samples" : 3,
              "libraries" : 4, "lanes"    : 5, "mapping" : 6,
              "trimming"  : 5}[show]

    for target in sorted(target for target in target_list if len(target) == length):
        if (show == "trimming") == (target[1] == "reads"):
            print ":".join(target)


def build_pipeline_targets(config, makefile):
    final_nodes = set()
    target_list = _make_target_list(config, makefile)
    for target in list(config.targets):
        key = tuple(target.split(":"))
        if key in target_list:
            final_nodes.update(target_list.get(key, ()))
            config.targets.remove(target)
    return final_nodes


def index_references(config, makefiles):
    references = {}
    for makefile in makefiles:
        for dd in makefile["Prefixes"].itervalues():
            reference = os.path.realpath(dd["Reference"])
            if reference not in references:
                references[reference] = \
                  MetaNode(description = "Reference Sequence",
                           dependencies = (FastaIndexNode(dd["Reference"]),
                                           BuildSequenceDictNode(config    = config,
                                                                 reference = dd["Reference"])))
            dd["Node"] = references[reference]


def parse_config(argv):
    config = ConfigParser.SafeConfigParser()
    config_paths = (os.path.join(os.path.expanduser('~'), ".pypeline.conf"),
                    "/etc/pypeline.conf")

    for config_path in config_paths:
        if os.path.exists(config_path):
            config.read(config_path)
            break

    try:
        defaults = dict(config.items("Defaults"))
    except ConfigParser.NoSectionError:
        defaults = {}

    parser = optparse.OptionParser()
    parser.add_option("--verbose", action = "store_true", default = defaults.get("verbose", False),
                      help = "Print the full dependency-tree every time a node is updated.")
    parser.add_option("--allow-missing-input-files", action = "store_true", default = False,
                      help = "Allow processing of lanes, even if the original input files are no-longer " \
                             "accesible, if for example a network drive is down. This option should be " \
                             "used with care!")

    group  = optparse.OptionGroup(parser, "Scheduling")
    group.add_option("--bowtie2-max-threads", type = int, default = defaults.get("bowtie2_max_threads", 4),
                     help = "Maximum number of threads to use per BWA instance [%default]")
    group.add_option("--bwa-max-threads", type = int, default = defaults.get("bwa_max_threads", 4),
                     help = "Maximum number of threads to use per BWA instance [%default]")
    group.add_option("--max-threads", type = int, default = defaults.get("max_threads", 14),
                     help = "Maximum number of threads to use in total [%default]")
    group.add_option("--dry-run", action = "store_true", default = False,
                     help = "If passed, only a dry-run in performed, the dependency tree is printed, and no tasks are executed.")
    parser.add_option_group(group)

    group  = optparse.OptionGroup(parser, "Required paths")
    group.add_option("--jar-root", default = os.path.expanduser(defaults.get("jar_root", os.path.join('~', "install", "jar_root"))),
                     help = "Folder containing Picard JARs (http://picard.sf.net), " \
                            "and GATK (www.broadinstitute.org/gatk). " \
                            "The latter is only required if realigning is enabled. " \
                            "[%default]")
    group.add_option("--temp-root", default = os.path.expanduser(defaults.get("temp_root", os.path.join('~', "scratch", "bam_pypeline"))),
                     help = "Location for temporary files and folders [%default/]")
    group.add_option("--destination", default = None,
                     help = "The destination folder for result files. By default, files will be "
                            "placed in the same folder as the makefile which generated it.")
    parser.add_option_group(group)

    group  = optparse.OptionGroup(parser, "Output files and orphan files")
    group.add_option("--list-output-files", action = "store_true", default = False,
                     help = "List all files generated by pipeline for the makefile(s).")
    group.add_option("--list-orphan-files", action = "store_true", default = False,
                     help = "List all files at destination not generated by the pipeline. " \
                            "Useful for cleaning up after making changes to a makefile.")
    parser.add_option_group(group)

    group  = optparse.OptionGroup(parser, "Targets")
    group.add_option("--target", dest = "targets", action = "append", default = [],
                     help = "Only execute nodes required to build specified target.")
    group.add_option("--list-targets", default = None,
                     help = "List all targets at a given resolution (target, sample, library, lane, reads)")
    parser.add_option_group(group)

    config, args = parser.parse_args(argv)

    config.targets = set(config.targets)
    targets_by_name = ("targets", "prefixes", "samples", "libraries", "lanes", "mapping", "trimming")
    if (config.list_targets is not None) and (config.list_targets not in targets_by_name):
        parser.error("ERROR: Invalid value for --list-targets (%s), valid values are '%s'." \
                     % (repr(config.list_targets), "', '".join(targets_by_name)))

    if config.list_output_files and config.list_orphan_files:
        parser.error("ERROR: Both --list-output-files and --list-orphan-files set!")
        return None

    return config, args


def walk_nodes(nodes, func, skip_nodes = None):
    if skip_nodes is None:
        skip_nodes = set()

    for node in nodes:
        if node in skip_nodes:
            continue
        elif not func(node):
            return False

        skip_nodes.add(node)
        if not walk_nodes(node.subnodes, func, skip_nodes):
            return False
        elif not walk_nodes(node.dependencies, func, skip_nodes):
            return False

    return True


def list_output_files(nodes):
    output_files = set()
    def collect_output_files(node):
        output_files.update(map(os.path.abspath, node.output_files))
        return True

    walk_nodes(nodes, collect_output_files)

    return output_files


def list_orphan_files(config, makefiles, nodes):
    files, mkfiles = set(), set()
    for mkfile in makefiles:
        mkfiles.add(os.path.abspath(mkfile["Statistics"]["Filename"]))
        for target in mkfile["Targets"]:
            glob_str = os.path.join(config.destination, target + "*")
            for root_filename in glob.glob(glob_str):
                if os.path.isdir(root_filename):
                    for (dirpath, _, filenames) in os.walk(root_filename):
                        files.update(os.path.abspath(os.path.join(dirpath, filename)) for filename in filenames)
                else:
                    files.add(os.path.abspath(root_filename))
    return (files - mkfiles) - list_output_files(nodes)


def main(argv):
    config_args = parse_config(argv)
    if not config_args:
        return 1

    config, args = config_args

    try:
        ui.print_info("Building BAM pipeline ...", file = sys.stderr)
        makefiles = read_makefiles(args)
        if not makefiles:
            ui.print_err("Plase specify at least one makefile!", file = sys.stderr)
            return 1
    except MakefileError, e:
        ui.print_err("Error reading makefile:\n\t%s" % \
                         "\n\t".join(str(e).split("\n")),
                         file = sys.stderr)
        return 1

    # Build .fai files for reference .fasta files
    index_references(config, makefiles)

    if config.list_targets:
        ui.print_info("Listing targets for %s ..." % config.list_targets, file = sys.stderr)
        for makefile in makefiles:
            list_targets_for(config, makefile, config.list_targets)
        return 0

    pipeline_func = build_pipeline_trimming
    if config.targets:
        pipeline_func = build_pipeline_targets
    elif os.path.basename(sys.argv[0]) != "trim_pipeline":
        pipeline_func = build_pipeline_full

    pipeline = pypeline.Pypeline(config)
    for makefile in makefiles:
        # If a destination is not specified, save results in same folder as makefile
        filename = makefile["Statistics"]["Filename"]
        old_destination = config.destination
        if old_destination is None:
            config.destination = os.path.dirname(filename)

        try:
            nodes = pipeline_func(config, makefile)
        except pypeline.node.NodeError, e:
            ui.print_err("Error while building pipeline for '%s':\n%s" % (filename, e))
            return 1

        config.destination = old_destination

        pipeline.add_nodes(nodes)

    if config.targets:
        ui.print_err("ERROR: Could not find --target(s): '%s'" % "', '".join(config.targets))
        ui.print_err("       Please use --list-targets to print list of valid target names.")
        return 1
    elif config.list_output_files:
        ui.print_info("Printing output files ...", file = sys.stderr)
        for filename in sorted(list_output_files(pipeline.nodes)):
            print(filename)
        return 0
    elif config.list_orphan_files:
        ui.print_info("Printing orphan files ...", file = sys.stderr)
        for filename in sorted(list_orphan_files(config, makefiles, pipeline.nodes)):
            print(filename)
        return 0

    ui.print_info("Running BAM pipeline ...", file = sys.stderr)
    if not pipeline.run(dry_run     = config.dry_run,
                        max_running = config.max_threads,
                        verbose     = config.verbose):
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
