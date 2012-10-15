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
import string
import optparse
import collections

import pypeline
import pypeline.ui as ui

from pypeline.nodes.bwa import SE_BWANode, PE_BWANode
from pypeline.nodes.gatk import IndelRealignerNode
from pypeline.nodes.picard import MergeSamFilesNode, MarkDuplicatesNode
from pypeline.nodes.coverage import CoverageNode
from pypeline.nodes.samtools import BAMIndexNode
from pypeline.nodes.adapterremoval import SE_AdapterRemovalNode, PE_AdapterRemovalNode

from pypeline.common.text import parse_padded_table
from pypeline.common.fileutils import swap_ext, add_postfix

import pypeline.tools.bam_pipeline as common
from pypeline.tools.bam_pipeline.nodes import *
from pypeline.tools.bam_pipeline.summary import SummaryTableNode
from pypeline.tools.bam_pipeline.validation import validate_records


_ADAPTERRM_SE_CACHE = {}
_ADAPTERRM_PE_CACHE = {}




def collect_records(filenames):
    records = collections.defaultdict(list)
    for filename in filenames:
        try:
            with open(filename) as mkfile:
                for record in parse_padded_table(mkfile):
                    record["files"] = common.paths.collect_files(record)
                    records[record["Name"]].append(record)
        except IOError, expt:
            ui.print_err("ERROR: Failed to read read makefile:\n\t%s\n" % (expt,))
            return None                        
                
    return records


def build_se_trim_nodes(config, bwa_prefix, record):
    try:
        trim_prefix, trim_node = _ADAPTERRM_SE_CACHE[id(record)]
    except KeyError:
        trim_prefix = os.path.join(common.paths.full_path(record, "reads"), "reads")
        trim_node   = SE_AdapterRemovalNode(input_files     = record["files"]["SE"], 
                                            output_prefix   = trim_prefix)
        _ADAPTERRM_SE_CACHE[id(record)] = (trim_prefix, trim_node)

    output_dir  = common.paths.full_path(record, bwa_prefix)
    output_file = os.path.join(output_dir, "reads.truncated.bam")
    
    read_group = "@RG\tID:%(Library)s\tSM:%(Sample)s\tLB:%(Library)s\tPU:%(Barcode)s\tPL:Illumina" % record
    aln  = SE_BWANode(input_file   = trim_prefix + ".truncated.gz",
                      output_file  = output_file,
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      min_quality  = config.bwa_min_quality,
                      threads      = config.bwa_max_threads,
                      dependencies = trim_node)

    return { "aligned"   : {"files" : [], 
                            "nodes" : []},
             "unaligned" : {"files" : [output_file],
                            "nodes" : [ValidateBAMFile(config, aln)]},
             "record"    : record}


def build_pe_trim_nodes(config, bwa_prefix, record):
    try:
        trim_prefix, trim_node = _ADAPTERRM_PE_CACHE[id(record)]
    except KeyError:
        trim_prefix = os.path.join(common.paths.full_path(record, "reads"), "reads")
        trim_node   = PE_AdapterRemovalNode(input_files_1 = record["files"]["PE_1"], 
                                            input_files_2 = record["files"]["PE_2"], 
                                            output_prefix = trim_prefix)
        _ADAPTERRM_PE_CACHE[id(record)] = (trim_prefix, trim_node)


    nodes = {}
    output_dir = common.paths.full_path(record, bwa_prefix)
    read_group = "@RG\tID:%(Library)s\tSM:%(Sample)s\tLB:%(Library)s\tPU:%(Barcode)s\tPL:Illumina" % record
    for filename in ("aln", "unaln"):
        input_filename  = trim_prefix + ".singleton.%s.truncated.gz" % (filename,)
        output_filename = os.path.join(output_dir, "reads.singleton.%s.truncated.bam" % (filename,))

        aln = SE_BWANode(input_file   = input_filename,
                         output_file  = output_filename,
                         prefix       = bwa_prefix,
                         read_group   = read_group,
                         min_quality  = config.bwa_min_quality,
                         threads      = config.bwa_max_threads,
                         dependencies = trim_node)
        nodes[filename] = [ValidateBAMFile(config, aln)]


    output_filename = os.path.join(output_dir, "reads.pairs.truncated.bam")
    aln  = PE_BWANode(input_file_1 = trim_prefix + ".pair1.truncated.gz",
                      input_file_2 = trim_prefix + ".pair2.truncated.gz",
                      output_file  = output_filename,
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      min_quality  = config.bwa_min_quality,
                      threads      = config.bwa_max_threads,
                      dependencies = trim_node)
    nodes["unaln"].append(ValidateBAMFile(config, aln))


    return { "aligned"   : {"files" : [os.path.join(output_dir, "reads.singleton.aln.truncated.bam")],
                            "nodes" : nodes["aln"]},
             "unaligned" : {"files" : [os.path.join(output_dir, "reads.singleton.unaln.truncated.bam"),
                                       os.path.join(output_dir, "reads.pairs.truncated.bam")],
                            "nodes" : nodes["unaln"]},
             "record"    : record }


def build_markdup_node(config, input_files, output_file, dependencies):
    mkdup_params = MarkDuplicatesNode.customize(config         = config, 
                                                input_bams     = input_files,
                                                output_bam     = output_file,
                                                output_metrics = swap_ext(output_file, ".metrics"),
                                                dependencies = dependencies)

    return ValidateBAMFile(config      = config,
                           node        = mkdup_params.build_node())


def build_kirdup_node(config, input_files, output_file, dependencies):
    node = FilterUniqueBAMNode(config       = config,
                               input_files  = input_files,
                               output_file  = output_file,
                               dependencies = dependencies)
    return ValidateBAMFile(config      = config,
                           node        = node) 


def build_library_merge_nodes(config, bwa_prefix, records):
    def_dict = lambda: {"aligned" : [], "unaligned" : []}
    nodes = collections.defaultdict(def_dict)
    files = collections.defaultdict(def_dict)
    paths = dict()
    for record in records:
        if common.paths.is_paired_end(record):
            entry = build_pe_trim_nodes(config, bwa_prefix, record)
        else:
            entry = build_se_trim_nodes(config, bwa_prefix, record)

        library = record["Library"]
        for key in ("aligned", "unaligned"):
            nodes[library][key].extend(entry[key]["nodes"])
            files[library][key].extend(entry[key]["files"])
        paths[record["Library"]] = common.paths.library_path(record, bwa_prefix)

    libraries = {}
    for library in nodes:
        if not files[library]["aligned"]:
            node = build_markdup_node(config        = config,
                                      input_files   = files[library]["unaligned"], 
                                      output_file   = paths[library] + ".unaligned.bam", 
                                      dependencies  = nodes[library]["unaligned"])
        else:
            markdup = build_markdup_node(config        = config,
                                         input_files   = files[library]["unaligned"], 
                                         output_file   = paths[library] + ".unaligned.markdup.bam", 
                                         dependencies  = nodes[library]["unaligned"])

            kirdup = build_kirdup_node(config, 
                                       input_files  = files[library]["aligned"], 
                                       output_file  = paths[library] + ".unaligned.kirdup.bam", 
                                       dependencies = nodes[library]["aligned"])

            merged = MergeSamFilesNode(config       = config,
                                       input_bams   = [paths[library] + ".unaligned.markdup.bam",
                                                       paths[library] + ".unaligned.kirdup.bam"],
                                       output_bam   = paths[library]  + ".unaligned.bam",
                                       dependencies = [markdup, kirdup])

            node = ValidateBAMFile(config      = config,
                                   node        = merged) 

            
        yield library, paths[library], node

      

def build_library_nodes(config, bwa_prefix, name, records):
    library_nodes = dict()
    
    for (library, library_prefix, nodes) in build_library_merge_nodes(config, bwa_prefix, records):
        output_file    = library_prefix + ".unaligned.bam"      
        nodes = MapDamageNode(reference        = common.paths.reference_sequence(bwa_prefix),
                              input_file       = output_file,
                              output_directory = common.paths.mapdamage_path(name, library, bwa_prefix),
                              dependencies     = nodes)

        library_nodes[library_prefix] = nodes
    return library_nodes



def build_merged_nodes(config, bwa_prefix, name, records):
    library_nodes = build_library_nodes(config, bwa_prefix, name, records)
    nodes = library_nodes.values()

    input_files = [(prefix + ".unaligned.bam") for prefix in library_nodes]
    output_file = common.paths.target_path(records[0], bwa_prefix)
    nodes = MergeSamFilesNode(config       = config, 
                              input_bams   = input_files,
                              output_bam   = output_file,
                              dependencies = nodes)

    log_file = common.paths.prefix_path(records[0], bwa_prefix) + ".unaligned.validated"
    nodes = ValidateBAMFile(config      = config,
                            node        = nodes,
                            log_file    = log_file)

    nodes = CoverageNode(input_file   = output_file,
                         name         = name,
                         dependencies = nodes)

    if not config.gatk_jar:
        return nodes

    unaligned_bam  = output_file
    output_file    = add_postfix(unaligned_bam, ".realigned")
    
    prefix         = common.paths.prefix_path(records[0], bwa_prefix)
    intervals_file = prefix + ".unaligned.intervals"
    validated_file = prefix + ".realigned.validated"
        
    nodes = IndelRealignerNode(config       = config, 
                               reference    = common.paths.reference_sequence(bwa_prefix),
                               infile       = unaligned_bam,
                               outfile      = output_file,
                               intervals    = intervals_file,
                               dependencies = nodes)

    nodes = ValidateBAMFile(config   = config,
                            node     = nodes,
                            log_file = validated_file)
    
    return CoverageNode(input_file   = output_file,
                    name         = name,
                    dependencies = nodes)


def build_nodes(config, bwa_prefix, name, records):
    targets = collections.defaultdict(list)
    any_pe_lanes = False
    for record in records:
        if not all(record["files"].values()):
            ui.print_err("Could not find files for record:")
            ui.print_err("\t- Name: %(Name)s\n\t- Sample: %(Sample)s\n\t- Library: %(Library)s\n\t- Barcode: %(Barcode)s" % record)
            if common.paths.is_paired_end(record):
                ui.print_err("\t- Found %i mate 1 files" % len(record["files"]["PE_1"]))
                ui.print_err("\t- Found %i mate 2 files" % len(record["files"]["PE_2"]))
                ui.print_err("Maybe this is a SE lane?")

            return None

    return build_merged_nodes(config, bwa_prefix, name, records)




def parse_config(argv):
    parser = optparse.OptionParser()
    parser.add_option("--destination", default = "./results",
                      help = "The destination folder for result files [%default/]")
    parser.add_option("--picard-root", default = None,
                      help = "Folder containing Picard JARs (http://picard.sf.net)")
    parser.add_option("--temp-root", default = "/tmp",
                      help = "Location for temporary files and folders [%default/]")
    parser.add_option("--bwa-prefix", action = "append", default = [],
                      help = "BWA prefix to align the input files against.")
    parser.add_option("--bwa-prefix-mito", default = None,
                      help = "BWA prefix of mitochondrial genome to align input files against (used in summary).")
    parser.add_option("--bwa-prefix-nuclear", default = None,
                      help = "BWA prefix of nuclear genome to align input files against (used in summary).")
    parser.add_option("--bwa-min-quality", default = 25, type = int,
                      help = "Minimum mapping quality (Phred score) of hits produced by BWA. " \
                             "Hits with a mapping quality below this value are filtered. [%default]")
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

    config.bwa_prefix = set(config.bwa_prefix)
    if config.bwa_prefix_mito:
        config.bwa_prefix.add(config.bwa_prefix_mito)
    if config.bwa_prefix_nuclear:
        config.bwa_prefix.add(config.bwa_prefix_nuclear)

    for prefix in config.bwa_prefix:
        if not common.paths.reference_sequence(prefix):
            errors.append("ERROR: Could not find reference sequence for prefix: '%s'" % prefix)
            errors.append("       Refernce sequences MUST have the extensions .fasta or .fa, and")
            errors.append("       be located at '${prefix}', '${prefix}.fa' or '${prefix}.fasta'.")
        elif (set(prefix) & set(string.whitespace)):
            errors.append("ERROR: BWA prefix must not contain whitespace:\n\t- Prefix: %s" % prefix)

        label = common.paths.prefix_to_filename(prefix)
        if label in ("*", "mito", "nuclear", "reads"):
            errors.append("ERROR: Prefix name is reserved keyword ('*', 'mito', 'nuclear', 'reads'), please rename:\n\t- Prefix: %s" % prefix)
    
    if not config.bwa_prefix:
        errors.append("ERROR: At least one BWA prefix must be specified using --bwa-prefix,")
        errors.append("       --bwa-prefix-mito, or --bwa-prefix-nucl.")

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
    common.paths.ROOT = config.destination

    records = collect_records(args)
    if records is None:
        return 1

    if not validate_records(records):
        return 1

    pipeline = pypeline.Pypeline(config)
    for (target, runs) in records.iteritems():
        nodes = []
        for bwa_prefix in config.bwa_prefix:
            current_nodes = build_nodes(config, bwa_prefix, target, runs)
            if not current_nodes:
                return 1
            nodes.append(current_nodes)

        pipeline.add_nodes(SummaryTableNode(config       = config,
                                            records      = runs,
                                            dependencies = nodes))

    if not pipeline.run(dry_run = config.dry_run, max_running = config.max_threads):
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

