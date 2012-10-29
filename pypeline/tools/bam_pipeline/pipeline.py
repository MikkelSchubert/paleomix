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
from pypeline.nodes.samtools import BAMIndexNode
from pypeline.nodes.adapterremoval import SE_AdapterRemovalNode, PE_AdapterRemovalNode

from pypeline.common.text import parse_padded_table
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
        node  = SE_AdapterRemovalNode(input_files   = reads["Raw"]["SE"],
                                      output_prefix = output_prefix)

        reads["Trimmed"] = {"Single" : output_prefix + ".truncated.gz"}
    else:
        node  = PE_AdapterRemovalNode(input_files_1 = reads["Raw"]["PE_1"],
                                      input_files_2 = reads["Raw"]["PE_2"],
                                      output_prefix = output_prefix)
        reads["Trimmed"] = {"Single"    : output_prefix + ".singleton.unaln.truncated.gz",
                            "Paired"    : output_prefix + ".pair{Pair}.truncated.gz",
                            "Collapsed" : output_prefix + ".singleton.aln.truncated.gz" }
    return node


def build_bwa_nodes(config, target, sample, library, barcode, record, dependencies):
    # Library is used as ID, to allow at-a-glance identification of
    # the source library for any given read in a BAM file.
    read_group = "@RG\tID:%s\tSM:%s\tLB:%s\tPU:%s\tPL:%s" \
        % (library, sample, library, barcode, record["Options"]["Platform"])

    reads = record["Reads"]
    reads["BAM"] = {}
    for (genome, prefix) in record["Prefixes"].iteritems():
        reads["BAM"][genome] = {}
        output_dir = os.path.join(config.destination, target, genome, sample, library, barcode)
        for (key, input_filename) in reads["Trimmed"].iteritems():
            # TODO: Make seeding optional
            options = ["minQ%i" % record["Options"]["BWA_MinQuality"]]
            if not record["Options"]["BWA_UseSeed"]:
                options.append("noSeed")

            output_filename = os.path.join(output_dir, "%s.%s.bam" \
                                               % (key.lower(), ".".join(options)))
            
            # Common BWA parameters
            parameters = {"output_file"  : output_filename,
                          "prefix"       : prefix["Path"],
                          "reference"    : prefix["Reference"],
                          "threads"      : config.bwa_max_threads,
                          "dependencies" : dependencies}

            if paths.is_paired_end(input_filename):
                params = PE_BWANode.customize(input_file_1 = input_filename.format(Pair = 1),
                                              input_file_2 = input_filename.format(Pair = 2),
                                              **parameters)
                params.commands["sampe"].set_parameter("-r", read_group)
                if not record["Options"]["BWA_UseSeed"]:
                    params.commands["aln_1"].set_parameter("-l", 2**16 - 1)
                    params.commands["aln_2"].set_parameter("-l", 2**16 - 1)
            else:
                params = SE_BWANode.customize(input_file   = input_filename,
                                              **parameters)
                params.commands["samse"].set_parameter("-r", read_group)
                if not record["Options"]["BWA_UseSeed"]:
                    params.commands["aln"].set_parameter("-l", 2**16 - 1)

            params.commands["filter"].set_parameter('-q', record["Options"]["BWA_MinQuality"])

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
    record = copy.deepcopy(record)
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
                          "LibCoverage"  : [lib_coverage]}

    return merged


def build_sample_nodes(config, target, sample, libraries):
    collected = collections.defaultdict(dict)
    for (library, barcodes) in libraries.iteritems():
        for (genome, record) in build_library_nodes(config, target, sample, library, barcodes).iteritems():
            collected[genome][library] = record

    return collected


def build_mapdamage_nodes(config, prefixes, target, genome, libraries):
    for (library, record) in libraries.iteritems():
        reference = prefixes[genome]["Reference"]
        output_directory = os.path.join(config.destination, "%s.%s.mapDamage" % (target, genome), library)
        yield MapDamageNode(reference        = reference,
                            input_file       = record["Filename"],
                            output_directory = output_directory,
                            dependencies     = record["Node"])


def build_target_nodes(config, prefixes, target, samples):
    input_records = collections.defaultdict(list)
    mapdamage_records = collections.defaultdict(list)
    for (sample, libraries) in samples.iteritems():
        for (genome, libraries) in build_sample_nodes(config, target, sample, libraries).iteritems():
            input_records[genome].extend(libraries.values())
            mapdamage_records[genome].extend(build_mapdamage_nodes(config, prefixes, target, genome, libraries))

    nodes = []
    for (genome, records) in input_records.iteritems():
        output_filename = os.path.join(config.destination, "%s.%s.bam" % (target, genome))
        node = MergeSamFilesNode(config       = config,
                                 input_bams   = [record["Filename"] for record in records],
                                 output_bam   = output_filename,
                                 dependencies = [record["Node"] for record in records])
        node = ValidateBAMFile(config      = config,
                               node        = node)

        small_coverage  = MetaNode(description = "Lanes and Libraries",
                                   subnodes    = sum((record["LaneCoverage"] + record["LibCoverage"] for record in records), []))



        
        library_coverage = sum((record["LibCoverage"] for record in records), [])
        large_coverage  = [MergeCoverageNode(input_files  = sum((list(node.output_files) for node in library_coverage), []),
                                             output_file  = swap_ext(output_filename, ".coverage"),
                                             dependencies = small_coverage)]

        if config.gatk_jar:
            aligned = IndelRealignerNode(config       = config,
                                         reference    = prefixes[genome]["Reference"],
                                         infile       = output_filename,
                                         outfile      = add_postfix(output_filename, ".realigned"),
                                         intervals    = os.path.join(config.destination, target, genome + ".intervals"),
                                         dependencies = node)
            node = ValidateBAMFile(config      = config,
                                   node        = aligned)
            large_coverage += [CoverageNode(input_file   = add_postfix(output_filename, ".realigned"),
                                            name         = target,
                                            dependencies = node)]

        coverage = MetaNode(description  = "Coverage",
                            dependencies = (small_coverage,
                                            MetaNode(description = "Final BAMs",
                                                     subnodes    = large_coverage)))
#        summary  = SummaryTableNode(config       = config,
#                                    prefixes     = prefixes,
#                                    target       = target,
#                                    samples      = samples,
#                                    dependencies = coverage)
        mapdamage = MetaNode(description = "MapDamage",
                             subnodes    = mapdamage_records[genome])
        statistics =  MetaNode(description  = "Statistics:",
                               dependencies = (coverage, mapdamage))


        nodes.append(MetaNode(description  = "Genome: %s" % genome,
                              dependencies = (node, statistics)))

    return MetaNode(description  = "Target: %s" % target,
                    dependencies = nodes)


def build_nodes(config, makefile):
    nodes = []
    for (target, samples) in makefile["Targets"].iteritems():
        nodes.append(build_target_nodes(config, makefile["Prefixes"], target, samples))
    return nodes


def parse_config(argv):
    parser = optparse.OptionParser()
    parser.add_option("--destination", default = "./results",
                      help = "The destination folder for result files [%default/]")
    parser.add_option("--picard-root", default = None,
                      help = "Folder containing Picard JARs (http://picard.sf.net)")
    parser.add_option("--temp-root", default = "/tmp",
                      help = "Location for temporary files and folders [%default/]")
    parser.add_option("--bwa-max-threads", type = int, default = 4,
                      help = "Maximum number of threads to use per BWA instance [%default]")
    parser.add_option("--max-threads", type = int, default = 14,
                      help = "Maximum number of threads to use in total [%default]")
    parser.add_option("--gatk-jar", default = None,
                      help = "Location of GenomeAnalysisTK.jar (www.broadinstitute.org/gatk)." \
                             "If specified, BAM files are realigned using the IndelRealigner tool.")
    parser.add_option("--dry-run", action = "store_true", default = False,
                      help = "If passed, only a dry-run in performed, and no tasks are executed.")
    config, args = parser.parse_args(argv)

    errors, warnings = [], []
    if not config.picard_root:
        errors.append("ERROR: --picard-root must be set to the location of the Picard JARs: This is required")
        errors.append("       for ValidateSamFile.jar, MarkDuplicates.jar, and possibly more.")
    elif not os.path.isdir(config.picard_root):
        errors.append("ERROR: Path passed to --picard_root is not a directory: %s" % config.picard_root)

    if not config.gatk_jar:
        warnings.append("WARNING: --gatk-jar not set, indel realigned bams will not be produced.")
    elif not os.path.isfile(config.gatk_jar):
        errors.append("ERROR: Path passed to --gatk-jar is not a file: %s" % config.gatk_jar)


    map(ui.print_warn, warnings)
    if errors:
        map(ui.print_err, errors)
        ui.print_err("See --help for more information")
        return None

    return config, args


def main(argv):
    ui.print_info("Running BAM pipeline ...")
    config_args = parse_config(argv)
    if not config_args:
        return 1

    config, args = config_args
    paths.ROOT = config.destination

    try:    
        makefiles = validate_makefiles(read_makefiles(args))
        if not makefiles:
            ui.print_err("Plase specify at least one makefile!")
            return 1
    except MakefileError, e:
        ui.print_err("Error reading makefile:\n\t%s" % \
                         "\n\t".join(str(e).split("\n")))
        return 1

    pipeline = pypeline.Pypeline(config)
    for makefile in makefiles:
        pipeline.add_nodes(build_nodes(config, makefile))

    if not pipeline.run(dry_run = config.dry_run, max_running = config.max_threads):
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

