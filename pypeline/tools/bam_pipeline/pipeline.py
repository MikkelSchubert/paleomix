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
import glob
import time
import logging

import pypeline
import pypeline.yaml
import pypeline.logger

from pypeline.common.console import \
    print_err, \
    print_info

from pypeline.pipeline import \
    Pypeline
from pypeline.node import \
    MetaNode
from pypeline.nodes.picard import \
    BuildSequenceDictNode
from pypeline.nodes.samtools import \
    FastaIndexNode
from pypeline.nodes.bwa import \
    BWAIndexNode
from pypeline.nodes.bowtie2 import \
    Bowtie2IndexNode

from pypeline.tools.bam_pipeline.makefile import \
    MakefileError, \
    read_makefiles

import pypeline.tools.bam_pipeline.parts as parts
import pypeline.tools.bam_pipeline.config as bam_config
import pypeline.tools.bam_pipeline.mkfile as bam_mkfile



def _add_extra_nodes(config, makefile, targets):
    for target in targets:
        parts.add_statistics_nodes(config, makefile, target)

    return targets


def build_pipeline_trimming(config, makefile):
    """Builds only the nodes required to produce trimmed reads.
    This reduces the required complexity of the makefile to a minimum."""

    nodes = []
    for prefix in makefile["Prefixes"].itervalues():
        for (_, samples) in makefile["Targets"].iteritems():
            for (_, libraries) in samples.iteritems():
                for (_, barcodes) in libraries.iteritems():
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
        for (_, prefix) in makefile["Prefixes"].iteritems():
            samples = []
            for (sample_name, library_records) in sample_records.iteritems():
                libraries = []
                for (library_name, barcode_records) in library_records.iteritems():
                    lanes = []
                    for (barcode, record) in barcode_records.iteritems():
                        lanes.append(parts.Lane(config, prefix, record, barcode))

                    if any(lane.bams for lane in lanes):
                        libraries.append(parts.Library(config, target_name, prefix, lanes, library_name))

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
        if (show == "trimming") and (target[1] != "reads"):
            continue
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
    references         = {}
    references_bwa     = {}
    references_bowtie2 = {}
    for makefile in makefiles:
        for dd in makefile["Prefixes"].itervalues():
            reference  = dd["Reference"]
            if reference not in references:
                faidx_node   = FastaIndexNode(dd["Reference"])
                dict_node    = BuildSequenceDictNode(config    = config,
                                                     reference = reference)
                bwa_node     = BWAIndexNode(input_file = reference)
                bowtie2_node = Bowtie2IndexNode(input_file = reference)

                references[reference] = \
                  MetaNode(description = "Reference Sequence",
                           dependencies = (faidx_node, dict_node))
                references_bwa[reference] = \
                  MetaNode(description = "Reference Sequence",
                           dependencies = (faidx_node, dict_node, bwa_node))
                references_bowtie2[reference] = \
                  MetaNode(description = "Reference Sequence",
                           dependencies = (faidx_node, dict_node, bowtie2_node))

            dd["Node"]         = references[reference]
            dd["Node:BWA"]     = references_bwa[reference]
            dd["Node:Bowtie2"] = references_bowtie2[reference]


def list_orphan_files(config, makefiles, pipeline):
    files, mkfiles = set(), set()
    for mkfile in makefiles:
        mkfile_path = mkfile["Statistics"]["Filename"]
        mkfiles.add(os.path.abspath(mkfile_path))
        for target in mkfile["Targets"]:
            destination = config.destination
            if not destination:
                destination = os.path.dirname(mkfile_path)

            glob_str = os.path.join(destination, target + "*")
            for root_filename in glob.glob(glob_str):
                if os.path.isdir(root_filename):
                    for (dirpath, _, filenames) in os.walk(root_filename):
                        for filename in filenames:
                            fpath = os.path.join(dirpath, filename)
                            files.add(os.path.abspath(fpath))
                else:
                    files.add(os.path.abspath(root_filename))
    return (files - mkfiles) - pipeline.list_output_files()


def run(config, args):
    if not os.path.exists(config.temp_root):
        try:
            os.makedirs(config.temp_root)
        except OSError, error:
            print_err("ERROR: Could not create temp root:\n\t%s" % (error,))
            return 1

    if not os.access(config.temp_root, os.R_OK | os.W_OK | os.X_OK):
        print_err("ERROR: Insufficient permissions for temp root: '%s'"
                  % (config.temp_root,))
        return 1

    try:
        print_info("Building BAM pipeline ...", file=sys.stderr)
        makefiles = read_makefiles(config, args)
    except (MakefileError, pypeline.yaml.YAMLError, IOError), error:
        print_err("Error reading makefiles:",
                  "\n  %s:\n   " % (error.__class__.__name__,),
                  "\n    ".join(str(error).split("\n")),
                  file=sys.stderr)
        return 1

    logfile_template = time.strftime("bam_pipeline.%Y%m%d_%H%M%S_%%02i.log")
    pypeline.logger.initialize(config, logfile_template)
    logger = logging.getLogger(__name__)

    # Build .fai files for reference .fasta files
    index_references(config, makefiles)

    if config.list_targets:
        logger.info("Listing targets for %s ...", config.list_targets)
        for makefile in makefiles:
            # If a destination is not specified, save results in same folder as
            # the makefile
            filename = makefile["Statistics"]["Filename"]
            old_destination = config.destination
            if old_destination is None:
                config.destination = os.path.dirname(filename)

            list_targets_for(config, makefile, config.list_targets)
            config.destination = old_destination
        return 0

    pipeline_func = build_pipeline_trimming
    if config.targets:
        pipeline_func = build_pipeline_targets
    elif os.path.basename(sys.argv[0]) != "trim_pipeline":
        pipeline_func = build_pipeline_full

    pipeline = Pypeline(config)
    for makefile in makefiles:
        # If a destination is not specified, save results in same folder as the
        # makefile
        filename = makefile["Statistics"]["Filename"]
        old_destination = config.destination
        if old_destination is None:
            config.destination = os.path.dirname(filename)

        try:
            nodes = pipeline_func(config, makefile)
        except pypeline.node.NodeError, error:
            logger.error("Error while building pipeline for '%s':\n%s",
                         filename, error)
            return 1

        config.destination = old_destination

        pipeline.add_nodes(nodes)

    if config.targets:
        logger.error("ERROR: Could not find --target(s): '%s'", "', '".join(config.targets))
        logger.error("       Please use --list-targets to print list of valid target names.")
        return 1
    elif config.list_output_files:
        logger.info("Printing output files ...")
        pipeline.print_output_files()
        return 0
    elif config.list_orphan_files:
        logger.info("Printing orphan files ...")
        for filename in sorted(list_orphan_files(config, makefiles, pipeline)):
            print(filename)
        return 0
    elif config.list_executables:
        logger.info("Printing required executables ...")
        pipeline.print_required_executables()
        return 0

    logger.info("Running BAM pipeline ...")
    if not pipeline.run(dry_run=config.dry_run,
                        max_running=config.max_threads,
                        progress_ui=config.progress_ui):
        return 1

    return 0


def _print_usage():
    basename = os.path.basename(sys.argv[0])
    if basename == "paleomix":
        basename = "bam_pipeline"

    print_info("BAM Pipeline %s\n" % (pypeline.__version__,))
    print_info("Usage:")
    print_info("  -- %s help           -- Display this message" % basename)
    print_info("  -- %s makefile [...] -- Generate makefile from 'SampleSheet.csv' files." % basename)
    print_info("  -- %s dryrun [...]   -- Perform dry run of pipeline on provided makefiles." % basename)
    print_info("     %s                   Equivalent to 'bam_pipeline run --dry-run [...]'." % (" " * len(basename),))
    print_info("  -- %s run [...]      -- Run pipeline on provided makefiles." % basename)


def main(argv):
    try:
        config, args = bam_config.parse_config(argv)
        if args and args[0].startswith("dry"):
            config.dry_run = True
    except bam_config.ConfigError, error:
        print_err(error)
        return 1

    commands = ("makefile", "mkfile", "run", "dry_run", "dry-run", "dryrun")
    if (len(args) == 0) or (args[0] not in commands):
        _print_usage()
        return 1
    elif args[0] in ("mkfile", "makefile"):
        return bam_mkfile.main(args[1:])
    elif not args[1:]:
        _print_usage()
        print_err("\nPlease specify at least one makefile!")
        return 1

    return run(config, args[1:])
