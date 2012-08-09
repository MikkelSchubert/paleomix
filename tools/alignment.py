#!/usr/bin/python
import os
import sys
import glob
import optparse
import collections

import pypeline
import pypeline.ui as ui
import pypeline.common.fileutils as fileutils
from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.nodes.bwa import SE_BWANode, PE_BWANode
from pypeline.nodes.picard import ValidateBAMNode, MergeSamFilesNode, MarkDuplicatesNode
from pypeline.nodes.adapterremoval import SE_AdapterRemovalNode, PE_AdapterRemovalNode


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


def bwa_prefix_name(bwa_prefix):
    return os.path.splitext(os.path.basename(bwa_prefix))[0]


def collect_files(record):
    """ """
    files = {}
    for end in ("R1", "R2"):
        files[end] = list(sorted(glob.glob(record["Path"] % { "Pair" : end })))
    return files


def collect_records(filenames):
    def _split(line):
        return filter(None, line.strip().split("\t"))

    records = collections.defaultdict(list)
    for filename in filenames:
        with open(filename) as mkfile:
            header = _split(mkfile.readline())

            for line in mkfile:
                record = dict(zip(header, _split(line)))
                record["files"] = collect_files(record)
                records[record["Name"]].append(record)
                
    return records



_ADAPTERRM_CACHE = {}
def build_SE_nodes(config, bwa_prefix, record):
    try:
        trim_prefix, trim_node = _ADAPTERRM_CACHE[id(record)]
    except KeyError:
        trim_prefix = os.path.join("results", record["Name"], record["Library"], record["Barcode"], "reads")
        trim_node   = SE_AdapterRemovalNode(record["files"]["R1"], trim_prefix)
        _ADAPTERRM_CACHE[id(record)] = (trim_prefix, trim_node)


    output_dir  = os.path.join("results", record["Name"] + "." + bwa_prefix_name(bwa_prefix), record["Library"], record["Barcode"])
    output_file = os.path.join(output_dir, "reads.truncated.bam")
    
    read_group = "@RG\tID:%(Library)s\tSM:%(Sample)s\tLB:%(Library)s\tPU:%(Barcode)s\tPL:Illumina" % record
    aln  = SE_BWANode(input_file   = trim_prefix + ".truncated.gz",
                      output_file  = output_file,
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      threads      = config.max_bwa_threads,
                      dependencies = trim_node)

    return { "aligned"   : [],
             "unaligned" : [output_file],
             "nodes"     : [aln],
             "record"    : record}


def build_PE_nodes(config, bwa_prefix, record):
    try:
        trim_prefix, trim_node = _ADAPTERRM_CACHE[id(record)]
    except KeyError:
        trim_prefix = os.path.join("results", record["Name"], record["Library"], record["Barcode"], "reads")
        trim_node   = PE_AdapterRemovalNode(input_files_1 = record["files"]["R1"], 
                                            input_files_2 = record["files"]["R2"], 
                                            output_prefix = trim_prefix)
        _ADAPTERRM_CACHE[id(record)] = (trim_prefix, trim_node)


    nodes = []
    output_dir = os.path.join("results", record["Name"] + "." + bwa_prefix_name(bwa_prefix), record["Library"], record["Barcode"])
    read_group = "@RG\tID:%(Library)s\tSM:%(Sample)s\tLB:%(Library)s\tPU:%(Barcode)s\tPL:Illumina" % record
    for filename in ("aln", "unaln"):
        input_filename  = trim_prefix + ".singleton.%s.truncated.gz" % (filename,)
        output_filename = os.path.join(output_dir, "reads.singleton.%s.truncated.bam" % (filename,))

        aln = SE_BWANode(input_file   = input_filename,
                         output_file  = output_filename,
                         prefix       = bwa_prefix,
                         read_group   = read_group,
                         threads      = config.max_bwa_threads,
                         dependencies = trim_node)
        val = ValidateBAMNode(config   = config, 
                              bamfile  = output_filename,
                              ignore   = ["MATE_NOT_FOUND"],
                              dependencies = aln)
        nodes.append(val)


    output_filename = os.path.join(output_dir, "reads.pairs.truncated.bam")
    aln  = PE_BWANode(input_file_1 = trim_prefix + ".pair1.truncated.gz",
                      input_file_2 = trim_prefix + ".pair2.truncated.gz",
                      output_file  = output_filename,
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      threads      = config.max_bwa_threads,
                      dependencies = trim_node)
    val = ValidateBAMNode(config   = config, 
                          bamfile  = output_filename,
                          ignore   = ["MATE_NOT_FOUND"],
                          dependencies = aln)
    nodes.append(val)


    return { "aligned"   : [os.path.join(output_dir, "reads.singleton.aln.truncated.bam")],
             "unaligned" : [os.path.join(output_dir, "reads.singleton.unaln.truncated.bam"),
                            os.path.join(output_dir, "reads.pairs.truncated.bam")],
             "nodes"     : nodes,
             "record"    : record }


def build_lib_merge_nodes(config, bwa_prefix, nodes):
    for key in ("aligned", "unaligned"):
        files, deps = [], []
        for node in nodes:
            files.extend(node[key])
            deps.extend(node["nodes"])
        if not files:
            continue

        record = node["record"]
        outdir = os.path.join("results", record["Name"] + "." + bwa_prefix_name(bwa_prefix))
        target = os.path.join(outdir, "%s.%s.bam" % (record["Library"], key))
        
        merge = MergeSamFilesNode(config       = config, 
                                  input_files  = files, 
                                  output_file  = target,
                                  dependencies = deps)

        val   = ValidateBAMNode(config       = config, 
                                bamfile      = target,
                                ignore       = ["MATE_NOT_FOUND"],
                                dependencies = merge)       

        yield (key, target, val)

    
def build_dedupe_node(config, key, filename, node):
    if key == "aligned":
        target = fileutils.add_postfix(filename, ".kirdup",)
        dedup = FilterUniqueBAMNode(input_file   = filename,
                                    output_file  = target,
                                    dependencies = node)
    elif key == "unaligned":
        target = fileutils.add_postfix(filename, ".markdup",)
        dedup = MarkDuplicatesNode(config       = config, 
                                   input_file   = filename,
                                   output_file  = target,
                                   dependencies = node)
    else:
        assert False


    val = ValidateBAMNode(config       = config, 
                          bamfile      = target,
                          ignore       = ["MATE_NOT_FOUND"],
                          dependencies = dedup)

    return filename, val

    


def build_nodes(config, bwa_prefix, sample, records):
    nodes = collections.defaultdict(list)
    for record in records:
        if record["files"]["R2"]:
            node = build_PE_nodes(config, bwa_prefix, record)
        elif record["files"]["R1"]:
            node = build_SE_nodes(config, bwa_prefix, record)
        else:
            assert False

        nodes[record["Library"]].append(node)
    
    deduped = []
    for nodes in nodes.values():
        for (key, filename, node) in build_lib_merge_nodes(config, bwa_prefix, nodes):
            deduped.append(build_dedupe_node(config, key, filename, node))


    target = os.path.join("results", "%s.%s.bam" % (sample, bwa_prefix_name(bwa_prefix)))
    return MergeSamFilesNode(config       = config,
                             input_files  = [filename for (filename, _) in deduped],
                             output_file  = target,
                             dependencies = [node for (_, node) in deduped])
            
    
    



def main(argv):
    parser = optparse.OptionParser()
    parser.add_option("--picard-root", None)
    parser.add_option("--temp-root", default = "/tmp")
    parser.add_option("--bwa-prefix", action = "append", default = [])
    parser.add_option("--max-threads", type = int, default = 14)
    parser.add_option("--max-bwa-threads", type = int, default = 4)
    parser.add_option("--run", action = "store_true", default = False)
    config, args = parser.parse_args(argv)
    pipeline = pypeline.Pypeline(config)

    if not config.picard_root:
        ui.print_err("--picard-root must be set to the location of the Picard JARs.")
        return 1

    records = collect_records(args)
    for bwa_prefix in config.bwa_prefix:
        for (sample, runs) in records.iteritems():
            pipeline.add_nodes(build_nodes(config, bwa_prefix, sample, runs))

    pipeline.run(dry_run = not config.run, max_running = config.max_threads)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

