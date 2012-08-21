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
import pypeline.common.fileutils as fileutils

from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.nodes.bwa import SE_BWANode, PE_BWANode
from pypeline.nodes.gatk import IndelRealignerNode
from pypeline.nodes.picard import ValidateBAMNode, MergeSamFilesNode, MarkDuplicatesNode
from pypeline.nodes.coverage import CoverageNode
from pypeline.nodes.samtools import BAMIndexNode
from pypeline.nodes.filesystem import HardlinkNode
from pypeline.nodes.adapterremoval import SE_AdapterRemovalNode, PE_AdapterRemovalNode

from pypeline.common.fileutils import swap_ext, add_postfix

import pypeline.tools.alignment_common as common
from pypeline.tools.alignment_common.summary import SummaryTableNode

_ADAPTERRM_SE_CACHE = {}
_ADAPTERRM_PE_CACHE = {}



class MapDamageNode(CommandNode):
    def __init__(self, reference, input_file, output_directory, dependencies):
        command = AtomicCmd(["mapDamage.pl", "map", "-c",
                             "-i", "%(IN_BAM)s",
                             "-d", "%(OUT_DIR)s",
                             "-r", reference],
                            IN_BAM  = input_file,
                            OUT_DIR = output_directory,
                            OUT_STDOUT = output_directory + ".stdout",
                            OUT_STDERR = output_directory + ".stderr")

        description =  "<mapDamage: '%s' -> '%s'>" % (input_file, output_directory)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


class FilterUniqueBAMNode(CommandNode):
    def __init__(self, input_file, output_file, dependencies = ()):
        command = AtomicCmd(["FilterUniqueBAM", "--PIPE", "--library"],
                            IN_STDIN   = input_file,
                            OUT_STDOUT = output_file)

        description =  "<FilterUniqueBAM: '%s' -> '%s'>" % (input_file, output_file)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)



def validate_records_unique(records):
    paths = collections.defaultdict(list)
    for runs in records.itervalues():
        for record in runs:
            paths[common.paths.full_path(record)].append(record)

    errors = False
    for records in paths.itervalues():
        if len(records) > 1:
            errors = True
            ui.print_err("ERROR: {0} too similar records found (combination of specified fields must be unique):".format(len(records)))
            ui.print_err("\t- Name:    {Name}\n\t\t- Sample:  {Sample}\n\t\t- Library: {Library}\n\t\t- Barcode: {Barcode}\n".format(**records[0]))

    return not errors


def validate_records_libraries(records):
    """Checks that library names are unique to each sample in a target,
    under the assumption that multiple libraries may be produced from
    a sample, but not vice versa."""
    errors = False
    for (name, runs) in sorted(records.items()):
        libraries = collections.defaultdict(set)
        for record in runs:
            libraries[record["Library"]].add(record["Sample"])

        for (library, samples) in sorted(libraries.items()):
            if len(samples) > 1:
                ui.print_err("ERROR: Multiple samples in one library:")
                ui.print_err("\t- Target:  {0}\n\t- Library: {1}\n\t- Samples: {2}\n".format(name, library, ", ".join(samples)))
                errors = True

    return not errors


def validate_records_paths(records):
    paths = collections.defaultdict(list)
    for runs in records.itervalues():
        for record in runs:
            template = record["Path"]
            current_paths = [template]
            if common.paths.is_paired_end(record):
                for end in (1, 2):
                    current_paths.append(template.format(Pair = end))
            
            for path in current_paths:
                paths[path].append(record)

    printed = list()
    for path in sorted(paths):
        current_records = tuple(sorted(paths[path]))
        if (len(current_records) > 1) and (current_records not in printed):
            printed.append(current_records)
            descriptions = []
            for (ii, record) in enumerate(current_records, start = 1):
                descriptions.append("\t- Record {0}:\n\t\t- Name:    {Name}\n\t\t- Sample:  {Sample}\n\t\t- Library: {Library}\n\t\t- Barcode: {Barcode}".format(ii, **record))
            ui.print_err("ERROR: Multiple records specify same path(s):\n{0}\n".format("\n".join(descriptions)))

    return not bool(printed)


def validate_records(records):
    """Tests for possible errors in the makefile, including
     - Non-unique entries (Name ... Barcode columns)
     - Paths entered multiple times
     - Libraries with multiple samples."""
    return validate_records_unique(records) \
        and validate_records_libraries(records) \
        and validate_records_paths(records)


def collect_records(filenames):
    def _split(line):
        return filter(None, line.strip().split())

    records = collections.defaultdict(list)
    for filename in filenames:
        with open(filename) as mkfile:
            header = None
            for line in mkfile:
                if line.startswith("#") or not line.strip():
                    continue
                elif not header:
                    header = _split(line)
                    continue
                
                record = dict(zip(header, _split(line)))
                record["files"] = common.paths.collect_files(record)
                records[record["Name"]].append(record)
                
    return records


def validate_bams(config, node, input_file = None, output_file = None):
    if not input_file:
        filenames = set()
        for filename in node.output_files:
            if filename.lower().endswith(".bai"):
                filenames.add(swap_ext(filename, ".bam"))
            elif filename.lower().endswith(".bam"):
                filenames.add(filename)
            
        assert len(filenames) == 1, filenames
        input_file = filenames.pop()

    return ValidateBAMNode(config       = config, 
                           bamfile      = input_file,
                           output_file  = output_file,
                           ignore       = ["MATE_NOT_FOUND"],
                           dependencies = node)


def merge_or_hardlink(config, input_files, output_file, dependencies, validation_log = None):
    if len(input_files) == 1:
        link_bam = HardlinkNode(input_file    = input_files[0],
                                output_file   = output_file,
                                dependencies  = dependencies)

        link_bai = HardlinkNode(input_file    = swap_ext(input_files[0], ".bai"),
                                output_file   = swap_ext(output_file,    ".bai"),
                                dependencies  = dependencies)
        
        return MetaNode(description  = str(link_bam), 
                        subnodes     = [link_bam, link_bai],
                        dependencies = dependencies)

    node = MergeSamFilesNode(config       = config, 
                             input_files  = input_files, 
                             output_file  = output_file,
                             dependencies = dependencies)
    
    return validate_bams(config      = config, 
                         node        = node,
                         output_file = validation_log)


def build_se_nodes(config, bwa_prefix, record):
    try:
        trim_prefix, trim_node = _ADAPTERRM_SE_CACHE[id(record)]
    except KeyError:
        trim_prefix = os.path.join(common.paths.full_path(record, "reads"), "reads")
        trim_node   = SE_AdapterRemovalNode(record["files"]["SE"], trim_prefix)
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
                            "nodes" : [validate_bams(config, aln)]},
             "record"    : record}


def build_pe_nodes(config, bwa_prefix, record):
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
        nodes[filename] = [validate_bams(config, aln)]


    output_filename = os.path.join(output_dir, "reads.pairs.truncated.bam")
    aln  = PE_BWANode(input_file_1 = trim_prefix + ".pair1.truncated.gz",
                      input_file_2 = trim_prefix + ".pair2.truncated.gz",
                      output_file  = output_filename,
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      min_quality  = config.bwa_min_quality,
                      threads      = config.bwa_max_threads,
                      dependencies = trim_node)
    nodes["unaln"].append(validate_bams(config, aln))


    return { "aligned"   : {"files" : [os.path.join(output_dir, "reads.singleton.aln.truncated.bam")],
                            "nodes" : nodes["aln"]},
             "unaligned" : {"files" : [os.path.join(output_dir, "reads.singleton.unaln.truncated.bam"),
                                       os.path.join(output_dir, "reads.pairs.truncated.bam")],
                            "nodes" : nodes["unaln"]},
             "record"    : record }


def build_lib_merge_nodes(config, bwa_prefix, nodes):
    for key in ("aligned", "unaligned"):
        files, deps = [], []
        for node in nodes:
            files.extend(node[key]["files"])
            deps.extend(node[key]["nodes"])
        if not files:
            continue

        record = nodes[0]["record"]
        outdir = common.paths.sample_path(record, bwa_prefix)
        target = os.path.join(outdir, "%s.%s.bam" % (record["Library"], key))
        
        merge  = MergeSamFilesNode(config       = config, 
                                   input_files  = files, 
                                   output_file  = target,
                                   dependencies = deps)
        validate = ValidateBAMNode(config       = config,
                                   bamfile      = target, 
                                   dependencies = merge)
        
        yield (key, target, validate)

    
def build_dedupe_node(config, key, filename, node):
    if key == "aligned":
        target = fileutils.add_postfix(filename, ".kirdup",)
        fltnode = FilterUniqueBAMNode(input_file   = filename,
                                   output_file  = target,
                                   dependencies = node)
        # Indexing of the file allows for rapid summation of hits when building the summary table.
        dedup   = BAMIndexNode(infile = target,
                               dependencies = fltnode)
    elif key == "unaligned":
        target = fileutils.add_postfix(filename, ".markdup",)
        dedup = MarkDuplicatesNode(config       = config, 
                                   input_file   = filename,
                                   output_file  = target,
                                   dependencies = node)
    else:
        assert False

    return target, validate_bams(config, dedup)

    


def build_nodes(config, bwa_prefix, name, records):
    nodes = collections.defaultdict(list)
    for record in records:
        if not all(record["files"].values()):
            ui.print_err("Could not find files for record:")
            ui.print_err("\t- Name: %(Name)s\n\t- Sample: %(Sample)s\n\t- Library: %(Library)s\n\t- Barcode: %(Barcode)s" % record)
            if common.paths.is_paired_end(record):
                ui.print_err("\t- Found %i mate 1 files" % len(record["files"]["PE_1"]))
                ui.print_err("\t- Found %i mate 2 files" % len(record["files"]["PE_2"]))
                ui.print_err("Maybe this is a SE lane?")

            return None
        elif common.paths.is_paired_end(record):
            node = build_pe_nodes(config, bwa_prefix, record)
        else:
            node = build_se_nodes(config, bwa_prefix, record)

        nodes[(record["Sample"], record["Library"])].append(node)
    
    deduped = collections.defaultdict(dict)
    for ((sample, library), nodes) in nodes.items():
        for (key, filename, node) in build_lib_merge_nodes(config, bwa_prefix, nodes):
            target, node = build_dedupe_node(config, key, filename, node)
            deduped[(sample, library)][target] = node

    merged = {}
    for ((sample, library), subnodes) in deduped.iteritems():
        record = { "Name"    : name, 
                   "Sample"  : sample,
                   "Library" : library }

        target = common.paths.library_path(record, bwa_prefix) + ".bam"
        assert target not in merged, target

        merge = merge_or_hardlink(config       = config,
                                  input_files  = subnodes.keys(),
                                  output_file  = target,
                                  dependencies = subnodes.values())

        damage = MapDamageNode(reference        = common.paths.reference_sequence(bwa_prefix),
                               input_file       = target,
                               output_directory = common.paths.mapdamage_path(record, bwa_prefix),
                               dependencies     = merge)
        
        merged[target] = damage


    record = { "Name"    : name, 
               "Sample"  : sample }
    target_aln   = common.paths.target_path(record, bwa_prefix)
    target_unaln = add_postfix(target_aln, ".unaligned")
    prefix = common.paths.prefix_to_filename(bwa_prefix)
    validated_target = os.path.join(config.destination, name, prefix + ".unaligned.validated")

    validated  = merge_or_hardlink(config         = config,
                               input_files    = merged.keys(),
                               output_file    = target_unaln,
                               dependencies   = merged.values(),
                               validation_log = validated_target)

    if config.gatk_jar:
        realigned_intervals = os.path.join(config.destination, name, prefix + ".unaligned.intervals")
        realigned = IndelRealignerNode(config       = config, 
                                       reference    = common.paths.reference_sequence(bwa_prefix),
                                       infile       = target_unaln,
                                       outfile      = target_aln,
                                       intervals    = realigned_intervals,
                                       dependencies = validated)
        validated_target = os.path.join(config.destination, name, prefix + ".aligned.validated")
        validated = validate_bams(config      = config, 
                                  node        = realigned, 
                                  input_file  = target_aln,
                                  output_file = validated_target)

        return CoverageNode(input_file   = target_aln,
                            name         = name,
                            dependencies = validated)
    else:
        return CoverageNode(input_file   = target_unaln,
                            name         = name,
                            dependencies = validated)
    



def parse_config(argv):
    parser = optparse.OptionParser()
    parser.add_option("--destination", default = "results")
    parser.add_option("--picard-root", None)
    parser.add_option("--temp-root", default = "/tmp")
    parser.add_option("--bwa-prefix", action = "append", default = [])
    parser.add_option("--bwa-prefix-mito", default = None)
    parser.add_option("--bwa-prefix-nuclear", default = None)
    parser.add_option("--bwa-min-quality", default = 25, type = int)
    parser.add_option("--bwa-max-threads", type = int, default = 4)
    parser.add_option("--max-threads", type = int, default = 14)
    parser.add_option("--gatk-jar")
    parser.add_option("--run", action = "store_true", default = False)
    config, args = parser.parse_args(argv)

    if not config.picard_root:
        ui.print_err("--picard-root must be set to the location of the Picard JARs.")
        return None
    elif not config.gatk_jar or not os.path.exists(config.gatk_jar):
        ui.print_warn("--gatk-jar not set, indel realigned bams will not be produced.")

    config.bwa_prefix = set(config.bwa_prefix)
    if config.bwa_prefix_mito:
        config.bwa_prefix.add(config.bwa_prefix_mito)
    if config.bwa_prefix_nuclear:
        config.bwa_prefix.add(config.bwa_prefix_nuclear)

    for prefix in config.bwa_prefix:
        if not common.paths.reference_sequence(prefix):
            ui.print_err("ERROR: Could not find reference sequence for prefix: '%s'" % prefix)
            ui.print_err("       Must have extension .fasta or .fa, and be located at '${prefix}', '${prefix}.fa' or '${prefix}.fasta'.")
            return None
        elif (set(prefix) & set(string.whitespace)):
            ui.print_err("ERROR: BWA prefix must not contain whitespace:\n\t- Prefix: %s" % prefix)
            return None

        label = common.paths.prefix_to_filename(prefix)
        if label in ("*", "mito", "nuclear", "reads"):
            ui.print_err("ERROR: Prefix name is reserved keyword ('*', 'mito', 'nuclear', 'reads'), please rename:\n\t- Prefix: %s" % prefix)
            return None
    
    return config, args

    


def main(argv):
    config_args = parse_config(argv)
    if not config_args:
        return 1

    config, args = config_args
    common.paths.ROOT = config.destination

    records = collect_records(args)
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
            

    pipeline.run(dry_run = not config.run, max_running = config.max_threads)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

