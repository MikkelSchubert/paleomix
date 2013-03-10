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

from pypeline.nodes.bwa import SE_BWANode, PE_BWANode
from pypeline.nodes.bowtie2 import Bowtie2Node
from pypeline.nodes.gatk import IndelRealignerNode
from pypeline.nodes.picard import BuildSequenceDictNode, \
                                  MergeSamFilesNode, \
                                  MarkDuplicatesNode
from pypeline.nodes.coverage import CoverageNode, MergeCoverageNode
from pypeline.nodes.samtools import FastaIndexNode, BAMIndexNode
from pypeline.nodes.adapterremoval import SE_AdapterRemovalNode, PE_AdapterRemovalNode

from pypeline.common.text import parse_padded_table
from pypeline.common.utilities import safe_coerce_to_tuple
from pypeline.common.fileutils import swap_ext, add_postfix, missing_files

import pypeline.tools.bam_pipeline.paths as paths
from pypeline.tools.bam_pipeline.nodes import *
from pypeline.tools.bam_pipeline.summary import SummaryTableNode
from pypeline.tools.bam_pipeline.makefile import *




def build_trimming_nodes(config, target, sample, library, barcode, record):
    output_prefix = os.path.join(config.destination, target, "reads", sample, library, barcode, "reads")

    reads = record["Reads"]

    if "SE" in reads["Raw"]:
        cmd = SE_AdapterRemovalNode.customize(input_files   = reads["Raw"]["SE"],
                                              output_prefix = output_prefix)

        reads["Trimmed"] = {"Single" : output_prefix + ".truncated.gz"}
    else:
        cmd = PE_AdapterRemovalNode.customize(input_files_1 = reads["Raw"]["PE_1"],
                                              input_files_2 = reads["Raw"]["PE_2"],
                                              output_prefix = output_prefix)

        reads["Trimmed"] = {"Single"    : output_prefix + ".singleton.unaln.truncated.gz",
                            "Paired"    : output_prefix + ".pair{Pair}.truncated.gz",
                            "Collapsed" : output_prefix + ".singleton.aln.truncated.gz" }

    if record["Options"]["QualityOffset"] == "Solexa":
        cmd.command.set_parameter("--qualitybase", 64)
    else:
        cmd.command.set_parameter("--qualitybase", record["Options"]["QualityOffset"])

    if any(missing_files(input_files) for input_files in reads["Raw"].itervalues()):
        if config.allow_missing_input_files:
            return ()

    return cmd.build_node()


def build_bwa_nodes(config, target, sample, library, barcode, record, dependencies):
    # Library is used as ID, to allow at-a-glance identification of
    # the source library for any given read in a BAM file.
    read_group = "@RG\tID:%s\tSM:%s\tLB:%s\tPU:%s\tPL:%s" \
        % (library, sample, library, barcode, record["Options"]["Platform"])

    reads = record["Reads"]
    reads["BAM"] = {}
    for (genome, prefix) in record["Prefixes"].iteritems():
        prefix_dependencies = []
        prefix_dependencies.extend(safe_coerce_to_tuple(dependencies))
        prefix_dependencies.extend(safe_coerce_to_tuple(prefix["Node"]))

        reads["BAM"][genome] = {}
        output_dir = os.path.join(config.destination, target, genome, sample, library, barcode)
        for (key, input_filename) in reads["Trimmed"].iteritems():
            options = ["minQ%i" % record["Options"]["Aligners"]["BWA"]["MinQuality"]]
            if not record["Options"]["Aligners"]["BWA"]["UseSeed"]:
                options.append("noSeed")

            output_filename = os.path.join(output_dir, "%s.%s.bam" \
                                               % (key.lower(), ".".join(options)))

            # Common BWA parameters
            parameters = {"output_file"  : output_filename,
                          "prefix"       : prefix["Path"],
                          "reference"    : prefix["Reference"],
                          "threads"      : config.bwa_max_threads,
                          "dependencies" : prefix_dependencies}

            if paths.is_paired_end(input_filename):
                params   = PE_BWANode.customize(input_file_1 = input_filename.format(Pair = 1),
                                                input_file_2 = input_filename.format(Pair = 2),
                                                **parameters)
                aln_keys, sam_key = ("aln_1", "aln_2"), "sampe"
            else:
                params   = SE_BWANode.customize(input_file   = input_filename,
                                              **parameters)
                aln_keys, sam_key = ("aln",), "samse"

            for aln_key in aln_keys:
                if not record["Options"]["Aligners"]["BWA"]["UseSeed"]:
                    params.commands[aln_key].set_parameter("-l", 2**16 - 1)
                if record["Options"]["QualityOffset"] in (64, "Solexa"):
                    params.commands[aln_key].set_parameter("-I")

            pg_tags  = "bwa:CL:%s" % " ".join(map(str, params.commands[aln_keys[0]].call)).replace("%", "%%")
            pg_tags += " | " + " ".join(map(str, params.commands[sam_key].call)).replace("%", "%%")
            params.commands["convert"].push_parameter("--update-pg-tag", pg_tags)

            params.commands[sam_key].set_parameter("-r", read_group)
            params.commands["filter"].set_parameter('-q', record["Options"]["Aligners"]["Bowtie2"]["MinQuality"])

            validate = ValidateBAMFile(config      = config,
                                       node        = params.build_node())
            coverage = CoverageNode(input_file   = output_filename,
                                    name         = target,
                                    dependencies = validate)

            reads["BAM"][genome][key] = {"Node"     : validate,
                                         "Filename" : output_filename,
                                         "LaneCoverage" : [coverage]}

    return record


def build_bowtie2_nodes(config, target, sample, library, barcode, record, dependencies):
    reads = record["Reads"]
    reads["BAM"] = {}
    for (genome, prefix) in record["Prefixes"].iteritems():
        prefix_dependencies = []
        prefix_dependencies.extend(safe_coerce_to_tuple(dependencies))
        prefix_dependencies.extend(safe_coerce_to_tuple(prefix["Node"]))

        reads["BAM"][genome] = {}
        output_dir = os.path.join(config.destination, target, genome, sample, library, barcode)
        for (key, input_filename) in reads["Trimmed"].iteritems():
            parameters = []
            options = ["minQ%i" % record["Options"]["Aligners"]["Bowtie2"]["MinQuality"]]
            output_filename = os.path.join(output_dir, "%s.%s.bam" \
                                               % (key.lower(), ".".join(options)))

            # Common BWA parameters
            parameters = {"output_file"  : output_filename,
                          "prefix"       : prefix["Path"],
                          "reference"    : prefix["Reference"],
                          "threads"      : config.bowtie2_max_threads,
                          "dependencies" : prefix_dependencies}

            if paths.is_paired_end(input_filename):
                params   = Bowtie2Node.customize(input_file_1 = input_filename.format(Pair = 1),
                                                 input_file_2 = input_filename.format(Pair = 2),
                                                 **parameters)
            else:
                params   = Bowtie2Node.customize(input_file_1 = input_filename,
                                                 input_file_2 = None,
                                                 **parameters)

            params.commands["filter"].set_parameter('-q', record["Options"]["Aligners"]["Bowtie2"]["MinQuality"])
            if record["Options"]["QualityOffset"] == 64:
                params.commands["aln"].set_parameter("--phred64")
            elif record["Options"]["QualityOffset"] == 33:
                params.commands["aln"].set_parameter("--phred33")
            else:
                params.commands["aln"].set_parameter("--solexa-quals")

            for (param, value) in record["Options"]["Aligners"]["Bowtie2"].iteritems():
                if (param != "MinQuality"):
                    for value in safe_coerce_to_tuple(value):
                        params.commands["aln"].push_parameter(param, value)

            pg_tags = "bowtie2:CL:%s" % " ".join(map(str, params.commands["aln"].call)).replace("%", "%%")
            params.commands["convert"].push_parameter("--update-pg-tag", pg_tags)

            # Library is used as ID, to allow at-a-glance identification of
            # the source library for any given read in a BAM file.
            params.commands["aln"].set_parameter("--rg-id", library)
            params.commands["aln"].push_parameter("--rg", "SM:%s" % sample)
            params.commands["aln"].push_parameter("--rg", "LB:%s" % library)
            params.commands["aln"].push_parameter("--rg", "PU:%s" % barcode)
            params.commands["aln"].push_parameter("--rg", "PL:%s" % record["Options"]["Platform"])
            params.commands["aln"].push_parameter("--rg", "PG:bowtie2")

            validate = ValidateBAMFile(config      = config,
                                       node        = params.build_node())
            coverage = CoverageNode(input_file   = output_filename,
                                    name         = target,
                                    dependencies = validate)

            reads["BAM"][genome][key] = {"Node"     : validate,
                                         "Filename" : output_filename,
                                         "LaneCoverage" : [coverage]}

    return record


def build_aln_nodes(config, target, sample, library, barcode, record, dependencies):
    if record["Options"]["Aligners"]["Program"] == "BWA":
        func = build_bwa_nodes
    elif record["Options"]["Aligners"]["Program"] == "Bowtie2":
        func = build_bowtie2_nodes
    else:
        assert "Aligner not implemented", record["Options"]["Aligners"]["Program"]

    return func(config, target, sample, library, barcode, record, dependencies)


def build_bam_cleanup_nodes(config, target, sample, library, barcode, record):
    tags = {"ID" : library, "SM" : sample, "LB" : library, "PU" : barcode,
            "PL" : record["Options"]["Platform"]}

    results = {}
    for (genome, alignments) in record["Reads"]["BAM"].iteritems():
        results[genome] = {}
        output_dir = os.path.join(config.destination, target, genome, sample, library, barcode)
        for (key, filename) in alignments.iteritems():
            output_filename = os.path.join(output_dir, "processed_%s.minQ%i.bam" \
                                               % (key.lower(), record["Options"]["Aligners"]["BWA"]["MinQuality"]))

            node = CleanupBAMNode(config    = config,
                                  reference = record["Prefixes"][genome]["Reference"],
                                  input_bam = filename,
                                  output_bam = output_filename,
                                  min_mapq   = record["Options"]["Aligners"]["BWA"]["MinQuality"],
                                  tags       = tags)
            node = ValidateBAMFile(config      = config,
                                   node        = node)

            if missing_files((filename,)) and config.allow_missing_input_files:
                node = ()

            coverage = CoverageNode(input_file   = output_filename,
                                    name         = target,
                                    dependencies = node)

            results[genome][key] = {"Node" : node,
                                    "Filename" : output_filename,
                                    "LaneCoverage" : [coverage]}

    record["Reads"]["BAM"] = results
    return record



def build_lane_nodes(config, target, sample, library, barcode, record):
    """

    """
    if "BAM" not in record["Reads"]:
        dependencies = ()
        if "Trimmed" not in record["Reads"]:
            dependencies = build_trimming_nodes(config, target, sample, library, barcode, record)

            # Allow specific read-types to be excluded from processsing
            for key in record["Options"]["ExcludeReads"]:
                record["Reads"]["Trimmed"].pop(key, None)

        return build_aln_nodes(config, target, sample, library, barcode, record, dependencies)

    return build_bam_cleanup_nodes(config, target, sample, library, barcode, record)


def build_rmduplicates_nodes(config, target, sample, library, input_records):
    results = {}
    for genome in input_records:
        results[genome] = {}

        collected_records = collections.defaultdict(list)
        for (key, records) in input_records[genome].iteritems():
            key = "kirdup" if (key == "Collapsed") else "markdup"
            collected_records[key].extend(records)

        for (key, cls) in (("kirdup", IndexedFilterUniqueBAMNode), ("markdup", MarkDuplicatesNode)):
            if key in collected_records:
                input_nodes = [record["Node"] for record in collected_records[key] if record["Node"]]
                input_files = [record["Filename"] for record in collected_records[key]]
                output_filename = os.path.join(config.destination, target, genome, sample, \
                                                   library + ".unaligned.%s.bam" % key)

                node = cls(config       = config,
                           input_bams   = input_files,
                           output_bam   = output_filename,
                           dependencies = input_nodes)
                node = ValidateBAMFile(config      = config,
                                       node        = node)

                cov  = CoverageNode(input_file   = output_filename,
                                    name         = target,
                                    dependencies = node)

                coverages = sum((record["LaneCoverage"] for record in collected_records[key]), [])
                results[genome][key] = [{"Node"     : node,
                                         "Filename" : output_filename,
                                         "LaneCoverage" : coverages,
                                         "LibCoverage"  : [cov]}]

    return results


def build_rescale_nodes(config, makefile, target, sample, library, input_records):
    results = {}
    for (genome, records_by_genome) in input_records.iteritems():
        results[genome] = {}

        records = []
        for records_by_key in records_by_genome.itervalues():
            records.extend(records_by_key)

        nodes       = [record.get("Node", []) for record in records]
        cov_lane    = sum((record.get("LaneCoverage", []) for record in records), [])
        cov_lib     = sum((record.get("LibCoverage", []) for record in records), [])
        input_files = [record.get("Filename", []) for record in records]
        reference   = makefile["Prefixes"][genome]["Reference"]

        output_filename = os.path.join(config.destination, target, genome, sample, \
                                       library + ".rescaled.bam")

        node = MapDamageRescaleNode(config       = config,
                                    reference    = reference,
                                    input_files  = input_files,
                                    output_file  = output_filename,
                                    dependencies = nodes)
        node = BAMIndexNode(infile         = output_filename,
                            dependencies   = node)
        node = ValidateBAMFile(config      = config,
                               node        = node)

        results[genome]["rescale"] = [{"Node"     : node,
                                       "Filename" : output_filename,
                                       "LaneCoverage" : cov_lib,
                                       "LibCoverage"  : cov_lib}]

    return results


def build_library_nodes(config, makefile, target, sample, library, barcodes):
    input_records = collections.defaultdict(lambda: collections.defaultdict(list))
    for (barcode, record) in barcodes.iteritems():
        record = build_lane_nodes(config, target, sample, library, barcode, record)

        for (genome, alignments) in record["Reads"]["BAM"].iteritems():
            for (key, alignment) in alignments.iteritems():
                input_records[genome][key].append(alignment)

    filter_duplicates = any(record["Options"]["PCRDuplicates"] for record in barcodes.itervalues())
    if filter_duplicates:
         input_records = build_rmduplicates_nodes(config, target, sample, library, input_records)

    if any(record["Options"]["RescaleQualities"] for record in barcodes.itervalues()):
        input_records = build_rescale_nodes(config, makefile, target, sample, library, input_records)

    merged = {}
    for genome in input_records:
        records = [record for records in input_records[genome].itervalues() for record in records]

        input_nodes = [record["Node"] for record in records]
        input_files = [record["Filename"] for record in records]
        input_cov_lane = sum((record["LaneCoverage"] for record in records), [])
        input_cov_lib  = input_cov_lane
        if filter_duplicates:
            input_cov_lib  = sum((record["LibCoverage"]  for record in records), [])

        node = MetaNode(description  = "Library: %s" % library,
                        dependencies = input_nodes)

        merged[genome] = {"Node"         : node,
                          "Filenames"    : input_files,
                          "LaneCoverage" : input_cov_lane,
                          "LibCoverage"  : input_cov_lib,
                          "Sample"       : sample,
                          "Library"      : library}

    return merged


def build_sample_nodes(config, makefile, target, sample, libraries):
    collected = collections.defaultdict(dict)
    for (library, barcodes) in libraries.iteritems():
        for (genome, record) in build_library_nodes(config, makefile, target, sample, library, barcodes).iteritems():
            collected[genome][library] = record

    return collected


def build_mapdamage_nodes(config, prefix, target, records):
    nodes = []
    for record in records:
        output_directory = os.path.join(config.destination, "%s.%s.mapDamage" % (target, prefix["Name"]), record["Library"])
        nodes.append(MapDamageNode(config           = config,
                                   reference        = prefix["Reference"],
                                   input_files      = record["Filenames"],
                                   output_directory = output_directory,
                                   dependencies     = record["Node"]))

    return MetaNode(description = "MapDamage",
                    subnodes    = nodes)


def build_statistics_nodes(config, makefile, prefix, target, records):
    output_filename = os.path.join(config.destination, "%s.%s.bam" % (target, prefix["Name"]))

    lane_coverage = sum((record["LaneCoverage"] for record in records), [])
    library_coverage = sum((record["LibCoverage"] for record in records), [])
    part_coverage = MetaNode(description = "Lanes and Libraries",
                              subnodes    = lane_coverage + library_coverage)
    full_coverage = [MergeCoverageNode(input_files  = sum((list(node.output_files) for node in library_coverage), []),
                                       output_file  = swap_ext(output_filename, ".coverage"),
                                       dependencies = part_coverage)]

    nodes = []

    if "Coverage" in makefile["Options"]["Features"]:
        nodes.append(MetaNode(description  = "Coverage",
                              dependencies = (part_coverage,
                                              MetaNode(description = "Final BAMs",
                                                       subnodes    = full_coverage))))
    elif "Summary" in makefile["Options"]["Features"]:
        nodes.append(MetaNode(description  = "Coverage",
                              dependencies = (part_coverage)))

    if "mapDamage" in makefile["Options"]["Features"]:
        nodes.append(build_mapdamage_nodes(config, prefix, target, records))

    if not nodes:
        return None

    return MetaNode(description  = "Statistics:",
                    dependencies = nodes)


def build_target_nodes(config, makefile, target, samples):
    prefixes = makefile["Prefixes"]
    input_records = collections.defaultdict(list)
    for (sample, libraries) in samples.iteritems():
        for (genome, libraries) in build_sample_nodes(config, makefile, target, sample, libraries).iteritems():
            input_records[genome].extend(libraries.values())

    nodes = []
    for (genome, records) in input_records.iteritems():
        library_files = sum((record["Filenames"] for record in records), [])
        library_nodes = [record["Node"] for record in records]

        library_meta  = MetaNode(description = "Libraries",
                                 dependencies = library_nodes)
        library_summary = MetaNode(description = "Libraries",
                                   subnodes    = library_nodes)

        final_nodes = []
        output_filename = os.path.join(config.destination, "%s.%s.bam" % (target, genome))

        if "Raw BAM" in makefile["Options"]["Features"]:
            node = MergeSamFilesNode(config       = config,
                                     input_bams   = library_files,
                                     output_bam   = output_filename,
                                     dependencies = library_summary)
            final_nodes.append(ValidateBAMFile(config      = config,
                                                node        = node))

        if ("Realigned BAM" in makefile["Options"]["Features"]):
            aligned = IndelRealignerNode(config       = config,
                                         reference    = prefixes[genome]["Reference"],
                                         infiles      = library_files,
                                         outfile      = add_postfix(output_filename, ".realigned"),
                                         intervals    = os.path.join(config.destination, target, genome + ".intervals"),
                                         dependencies = library_summary)
            final_nodes.append(ValidateBAMFile(config      = config,
                                               node        = aligned))


        meta_nodes = [library_meta]
        if final_nodes:
            meta_nodes.append(MetaNode(description  = "Final Nodes",
                                       dependencies = final_nodes))

        stats_meta = build_statistics_nodes(config, makefile, prefixes[genome], target, records)
        if stats_meta:
            meta_nodes.append(stats_meta)

        nodes.append(MetaNode(description  = "Genome: %s" % genome,
                              dependencies = meta_nodes))

    if "Summary" in makefile["Options"]["Features"]:
        nodes = SummaryTableNode(config       = config,
                                 prefixes     = prefixes,
                                 makefile     = makefile["Statistics"],
                                 target       = target,
                                 samples      = samples,
                                 records      = input_records,
                                 dependencies = nodes)

    return MetaNode(description  = "Target: %s" % target,
                    dependencies = nodes)


def build_pipeline_trimming(config, makefile):
    """Builds only the nodes required to produce trimmed reads.
    This reduces the required complexity of the makefile to a minimum."""

    nodes = []
    for (target, samples) in makefile["Targets"].iteritems():
        for (sample, libraries) in samples.iteritems():
            for (library, barcodes) in libraries.iteritems():
                for (barcode, record) in barcodes.iteritems():
                    if ("BAM" in record["Reads"]) or ("Trimmed" in record["Reads"]):
                        continue

                    nodes.append(build_trimming_nodes(config, target, sample, library, barcode, record))
    return nodes



def build_pipeline_full(config, makefile):
    nodes = []
    for (target, samples) in makefile["Targets"].iteritems():
        nodes.append(build_target_nodes(config, makefile, target, samples))
    return nodes


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
    group.add_option("--target", action = "append", default = [],
                     help = "Only execute nodes required to build specified files.")
    group.add_option("--list-output-files", action = "store_true", default = False,
                     help = "List all files generated by pipeline for the makefile(s).")
    group.add_option("--list-orphan-files", action = "store_true", default = False,
                     help = "List all files at destination not generated by the pipeline. " \
                            "Useful for cleaning up after making changes to a makefile.")
    parser.add_option_group(group)

    config, args = parser.parse_args(argv)

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


def collect_target_nodes(nodes, target_files):
    collected_nodes = set()
    target_files = set(map(os.path.realpath, target_files))
    def collect_nodes(node):
        for filename in node.output_files:
            filename = os.path.realpath(filename)
            if filename in target_files:
                target_files.remove(filename)
                collected_nodes.add(node)

        return bool(target_files)

    walk_nodes(nodes, collect_nodes)
    return collected_nodes


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

    pipeline_func = build_pipeline_trimming
    if os.path.basename(sys.argv[0]) != "trim_pipeline":
        pipeline_func = build_pipeline_full
        # Build .fai files for reference .fasta files
        index_references(config, makefiles)

    pipeline = pypeline.Pypeline(config)
    for makefile in makefiles:
        # If a destination is not specified, save results in same folder as makefile
        old_destination = config.destination
        if old_destination is None:
            config.destination = os.path.dirname(makefile["Statistics"]["Filename"])

        nodes = pipeline_func(config, makefile)
        config.destination = old_destination

        if config.target:
            nodes = collect_target_nodes(nodes, config.target)
        pipeline.add_nodes(nodes)

    if config.list_output_files:
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
