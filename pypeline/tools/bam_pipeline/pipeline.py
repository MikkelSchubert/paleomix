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

import pypeline
import pypeline.ui as ui

from pypeline.nodes.bwa import SE_BWANode, PE_BWANode
from pypeline.nodes.gatk import IndelRealignerNode
from pypeline.nodes.picard import MergeSamFilesNode, MarkDuplicatesNode
from pypeline.nodes.coverage import CoverageNode, MergeCoverageNode
from pypeline.nodes.samtools import FastaIndexNode
from pypeline.nodes.adapterremoval import SE_AdapterRemovalNode, PE_AdapterRemovalNode

from pypeline.common.text import parse_padded_table
from pypeline.common.utilities import safe_coerce_to_tuple
from pypeline.common.fileutils import swap_ext, add_postfix

import pypeline.tools.bam_pipeline.paths as paths
from pypeline.tools.bam_pipeline.nodes import *
from pypeline.tools.bam_pipeline.summary import SummaryTableNode
from pypeline.tools.bam_pipeline.makefile import *
from pypeline.tools.bam_pipeline.validation import validate_makefiles




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

    cmd.command.set_parameter("--qualitybase", record["Options"]["QualityOffset"])
    
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
            # TODO: Make seeding optional
            options = ["minQ%i" % record["Options"]["BWA_MinQuality"]]
            if not record["Options"]["BWA_UseSeed"]:
                options.append("noSeed")

            max_edit = record["Options"]["BWA_MaxEdit"]
            if record["Options"]["BWA_MaxEdit"] != 0.04:
                options.append("maxEdit%s" % max_edit)

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

            params.commands[sam_key].set_parameter("-r", read_group)
            params.commands["filter"].set_parameter('-q', record["Options"]["BWA_MinQuality"])
            
            for aln_key in aln_keys:
                params.commands[aln_key].set_parameter("-n", max_edit)
                if not record["Options"]["BWA_UseSeed"]:
                    params.commands[aln_key].set_parameter("-l", 2**16 - 1)
                if record["Options"]["QualityOffset"] == 64:
                    params.commands[aln_key].set_parameter("-I")

            validate = ValidateBAMFile(config      = config,
                                       node        = params.build_node())
            coverage = CoverageNode(input_file   = output_filename,
                                    name         = target,
                                    dependencies = validate)

            reads["BAM"][genome][key] = {"Node"     : validate, 
                                         "Filename" : output_filename,
                                         "LaneCoverage" : coverage}
    
    return record


def build_bam_cleanup_nodes(config, target, sample, library, barcode, record):
    tags = {"ID" : library, "SM" : sample, "LB" : library, "PU" : barcode,
            "PL" : record["Options"]["Platform"]}

    results = {}
    for (genome, alignments) in record["Reads"]["BAM"].iteritems():
        results[genome] = {}
        output_dir = os.path.join(config.destination, target, genome, sample, library, barcode)
        for (key, filename) in alignments.iteritems():
            output_filename = os.path.join(output_dir, "processed_%s.minQ%i.bam" \
                                               % (key.lower(), record["Options"]["BWA_MinQuality"]))

            node = CleanupBAMNode(config    = config,
                                  reference = record["Prefixes"][genome]["Reference"],
                                  input_bam = filename,
                                  output_bam = output_filename,
                                  min_mapq   = record["Options"]["BWA_MinQuality"],
                                  tags       = tags)
            node = ValidateBAMFile(config      = config,
                                   node        = node)
            coverage = CoverageNode(input_file   = output_filename,
                                    name         = target,
                                    dependencies = node)

            results[genome][key] = {"Node" : node,
                                    "Filename" : output_filename, 
                                    "LaneCoverage" : coverage}

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

        return build_bwa_nodes(config, target, sample, library, barcode, record, dependencies)

    return build_bam_cleanup_nodes(config, target, sample, library, barcode, record)


def build_rmduplicates_nodes(config, target, sample, library, input_records):
    results = {}
    for genome in input_records:
        results[genome] = {}

        collected_records = collections.defaultdict(list)
        for (key, records) in input_records[genome].iteritems():
            key = "kirdup" if (key == "Collapsed") else "markdup"
            collected_records[key].extend(records)

        for (key, cls) in (("kirdup", FilterUniqueBAMNode), ("markdup", MarkDuplicatesNode)):
            if key in collected_records:
                input_nodes = [record["Node"] for record in collected_records[key]]
                input_files = [record["Filename"] for record in collected_records[key]]
                output_filename = os.path.join(config.destination, target, genome, sample, \
                                                   library + ".unaligned.%s.bam" % key)

                node = cls(config       = config,
                           input_bams   = input_files,
                           output_bam   = output_filename,
                           dependencies = input_nodes)
                node = ValidateBAMFile(config      = config,
                                       node        = node)
                
                coverages = [record["LaneCoverage"] for record in collected_records[key]]
                results[genome][key] = [{"Node"     : node,
                                         "Filename" : output_filename,
                                         "LaneCoverage" : coverages}]

    return results


def build_library_nodes(config, target, sample, library, barcodes):
    filter_duplicates = False
    input_records = collections.defaultdict(lambda: collections.defaultdict(list))
    for (barcode, record) in barcodes.iteritems():
        record = build_lane_nodes(config, target, sample, library, barcode, record)
        filter_duplicates = record["Options"]["PCRDuplicates"]

        for (genome, alignments) in record["Reads"]["BAM"].iteritems():
            for (key, alignment) in alignments.iteritems():
                input_records[genome][key].append(alignment)

    if filter_duplicates:
        input_records = build_rmduplicates_nodes(config, target, sample, library, input_records)

    merged = {}
    for genome in input_records:
        input_nodes = [record["Node"] for records in input_records[genome].itervalues() for record in records]
        input_files = [record["Filename"] for records in input_records[genome].itervalues() for record in records]

        output_filename = os.path.join(config.destination, target, genome, sample, library + ".unaligned.bam")
        node = MergeSamFilesNode(config       = config,
                                 input_bams   = input_files,
                                 output_bam   = output_filename,
                                 dependencies = input_nodes)
        node = ValidateBAMFile(config      = config,
                               node        = node)

        lane_coverage = sum((record["LaneCoverage"] for records in input_records[genome].itervalues() for record in records), [])
        lib_coverage  = CoverageNode(input_file   = output_filename,
                                     name         = target,
                                     dependencies = node)

        node = MetaNode(description  = "Library: %s" % library,
                        dependencies = node)

        merged[genome] = {"Node"         : node,
                          "Filename"     : output_filename,
                          "LaneCoverage" : lane_coverage,
                          "LibCoverage"  : [lib_coverage],
                          "Sample"       : sample,
                          "Library"      : library}

    return merged


def build_sample_nodes(config, target, sample, libraries):
    collected = collections.defaultdict(dict)
    for (library, barcodes) in libraries.iteritems():
        for (genome, record) in build_library_nodes(config, target, sample, library, barcodes).iteritems():
            collected[genome][library] = record

    return collected


def build_mapdamage_nodes(config, prefix, target, records):
    nodes = []
    for record in records:
        output_directory = os.path.join(config.destination, "%s.%s.mapDamage" % (target, prefix["Name"]), record["Library"])
        nodes.append(MapDamageNode(reference        = prefix["Reference"],
                                   input_file       = record["Filename"],
                                   output_directory = output_directory,
                                   dependencies     = record["Node"]))

    return MetaNode(description = "MapDamage",
                    subnodes    = nodes)


def build_statistics_nodes(config, prefix, target, records, aligned_node):
    output_filename = os.path.join(config.destination, "%s.%s.bam" % (target, prefix["Name"]))

    lane_coverage = sum((record["LaneCoverage"] for record in records), [])
    library_coverage = sum((record["LibCoverage"] for record in records), []) 
    part_coverage = MetaNode(description = "Lanes and Libraries",
                              subnodes    = lane_coverage + library_coverage)
    full_coverage = [MergeCoverageNode(input_files  = sum((list(node.output_files) for node in library_coverage), []),
                                       output_file  = swap_ext(output_filename, ".coverage"),
                                       dependencies = part_coverage)]

    if aligned_node:
        full_coverage.append(CoverageNode(input_file   = add_postfix(output_filename, ".realigned"),
                                          name         = target,
                                          dependencies = aligned_node))

    coverage  = MetaNode(description  = "Coverage",
                         dependencies = (part_coverage,
                                         MetaNode(description = "Final BAMs",
                                                  subnodes    = full_coverage)))
    mapdamage = build_mapdamage_nodes(config, prefix, target, records)                  

    return MetaNode(description  = "Statistics:",
                    dependencies = (coverage, mapdamage))


def build_target_nodes(config, makefile, prefixes, target, samples):
    input_records = collections.defaultdict(list)
    for (sample, libraries) in samples.iteritems():
        for (genome, libraries) in build_sample_nodes(config, target, sample, libraries).iteritems():
            input_records[genome].extend(libraries.values())

    nodes = []
    for (genome, records) in input_records.iteritems():
        output_filename = os.path.join(config.destination, "%s.%s.bam" % (target, genome))
        node = MergeSamFilesNode(config       = config,
                                 input_bams   = [record["Filename"] for record in records],
                                 output_bam   = output_filename,
                                 dependencies = [record["Node"] for record in records])
        node = ValidateBAMFile(config      = config,
                               node        = node)

        aligned = None
        if config.gatk_jar:
            aligned = IndelRealignerNode(config       = config,
                                         reference    = prefixes[genome]["Reference"],
                                         infile       = output_filename,
                                         outfile      = add_postfix(output_filename, ".realigned"),
                                         intervals    = os.path.join(config.destination, target, genome + ".intervals"),
                                         dependencies = node)
            aligned = ValidateBAMFile(config      = config,
                                      node        = aligned)
            
        statistics = build_statistics_nodes(config, prefixes[genome], target, records, aligned)
        nodes.append(MetaNode(description  = "Genome: %s" % genome,
                              dependencies = (aligned or node, statistics)))

    summary  = SummaryTableNode(config       = config,
                                prefixes     = prefixes,
                                makefile     = makefile,
                                target       = target,
                                samples      = samples,
                                records      = input_records,
                                dependencies = nodes)

    return MetaNode(description  = "Target: %s" % target,
                    dependencies = summary)


def build_nodes(config, makefile):
    nodes = []
    for (target, samples) in makefile["Targets"].iteritems():
        nodes.append(build_target_nodes(config, makefile["Makefile"], makefile["Prefixes"], target, samples))
    return nodes


def index_references(makefiles):
    references = {}
    for makefile in makefiles:
        for dd in makefile["Prefixes"].itervalues():
            reference = os.path.realpath(dd["Reference"])
            if reference not in references:
                references[reference] = FastaIndexNode(dd["Reference"])
            dd["Node"] = references[reference]


def parse_config(argv):
    parser = optparse.OptionParser()
    parser.add_option("--target", action = "append", default = [],
                      help = "Only execute nodes required to build specified files.")
    parser.add_option("--list-output-files", action = "store_true", default = False,
                      help = "List all files generated by pipeline.")
    parser.add_option("--list-orphan-files", action = "store_true", default = False,
                      help = "List all files at destination not generated by the pipeline.")

    parser.add_option("--destination", default = "./results",
                      help = "The destination folder for result files [%default/]")
    parser.add_option("--picard-root", default = os.path.join(os.path.expanduser('~'), "install", "picard-tools"),
                      help = "Folder containing Picard JARs (http://picard.sf.net)")
    parser.add_option("--temp-root", default = "/tmp/" + pwd.getpwuid(os.getuid()).pw_name,
                      help = "Location for temporary files and folders [%default/]")
    parser.add_option("--bwa-max-threads", type = int, default = 4,
                      help = "Maximum number of threads to use per BWA instance [%default]")
    parser.add_option("--max-threads", type = int, default = 14,
                      help = "Maximum number of threads to use in total [%default]")
    parser.add_option("--gatk-jar", default = os.path.join(os.path.expanduser('~'), "install", "GATK", "GenomeAnalysisTK.jar"),
                      help = "Location of GenomeAnalysisTK.jar (www.broadinstitute.org/gatk)." \
                             "If specified, BAM files are realigned using the IndelRealigner tool.")
    parser.add_option("--dry-run", action = "store_true", default = False,
                      help = "If passed, only a dry-run in performed, and no tasks are executed.")
    parser.add_option("--non-verbose", action = "store_true", default = False,
                      help = "Only print running nodes while running (useful for large projects).")
    config, args = parser.parse_args(argv)

    errors, warnings = [], []
    if not os.path.exists(config.picard_root):
        errors.append("ERROR: Path passed to --picard_root does not exist: %s" % config.picard_root)
    elif not os.path.isfile(os.path.join(config.picard_root, "ValidateSamFile.jar")):
        errors.append("ERROR: Path passed to --picard-root does not appear to contain required JARs: %s" % config.picard_root)

    if not os.path.exists(config.gatk_jar):
        warnings.append("WARNING: GATK jar does not exist, indel realigned bams will not be produced.")
        config.gatk_jar = None
    elif not os.path.isfile(config.gatk_jar):
        errors.append("ERROR: Path passed to --gatk-jar is not a file: %s" % config.gatk_jar)

    if config.list_output_files and config.list_orphan_files:
        errors.append("ERROR: Both --list-output-files and --list-orphan-files set!")

    for warning in warnings:
        ui.print_warn(warning, file = sys.stderr)

    if errors:
        errors.append("See --help for more information")
        for error in errors:
            ui.print_err(warning, file = sys.stderr)
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
        output_files.update(node.output_files)
        return True

    walk_nodes(nodes, collect_output_files)

    return output_files


def list_orphan_files(config, makefiles, nodes):
    files = set()
    for mkfile in makefiles:
        for target in mkfile["Targets"]:
            for root_filename in glob.glob(os.path.join(config.destination, target + "*")):
                if os.path.isdir(root_filename):
                    for (dirpath, _, filenames) in os.walk(root_filename):
                        files.update(os.path.join(dirpath, filename) for filename in filenames)
    return files - list_output_files(nodes)


def main(argv):
    config_args = parse_config(argv)
    if not config_args:
        return 1

    config, args = config_args
    paths.ROOT = config.destination

    try:
        ui.print_info("Building BAM pipeline ...", file = sys.stderr)
        makefiles = validate_makefiles(read_makefiles(args))
        if not makefiles:
            ui.print_err("Plase specify at least one makefile!", file = sys.stderr)
            return 1
    except MakefileError, e:
        ui.print_err("Error reading makefile:\n\t%s" % \
                         "\n\t".join(str(e).split("\n")),
                         file = sys.stderr)
        return 1

    index_references(makefiles)
    pipeline = pypeline.Pypeline(config)
    for makefile in makefiles:
        nodes = build_nodes(config, makefile)
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
                        verbose     = not config.non_verbose):
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
