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


def collect_files(root, record):
    """ """
    files = {}
    for end in ("R1", "R2"):
        template = "%%(SampleID)s_%%(Index)s_L*%%(Lane)s_%s_*.fastq.gz" % (end,)
        files[end] = list(sorted(glob.glob(os.path.join(root, template % record))))
    return files


def read_alignment_records(root):
    """ """
    with open(os.path.join(root, "SampleSheet.csv")) as records:
        header = records.readline().strip().split(",")
        for line in records:
            record = dict(zip(header, line.strip().split(",")))
            record["files"] = collect_files(root, record)

            yield record


def collect_records(roots):        
    records = collections.defaultdict(list)
    for root in roots:
        for record in read_alignment_records(root):
            records[record["SampleID"]].append(record)
    return records


def build_SE_nodes(bwa_prefix, record):
    read_group = "@RG\tID:%(Index)s\tSM:%(SampleID)s\tLB:%(Index)s\tPU:%(FCID)s\tPL:Illumina" % record
    outdir = os.path.join("results", record["SampleID"] + "." + bwa_prefix_name(bwa_prefix), record["Index"], record["FCID"])
    prefix = os.path.join(outdir, "reads.truncated")

    trim = SE_AdapterRemovalNode(record["files"]["R1"], os.path.join(outdir, "reads"))
    aln  = SE_BWANode(input_file   = prefix + ".gz",
                      output_file  = prefix + ".bam",
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      dependencies = trim)

    return { "aligned"   : [],
             "unaligned" : [prefix + ".bam"],
             "nodes"     : [aln],
             "record"    : record}


def build_PE_nodes(config, bwa_prefix, record):
    read_group = "@RG\tID:%(Index)s\tSM:%(SampleID)s\tLB:%(Index)s\tPU:%(FCID)s\tPL:Illumina" % record
    outdir = os.path.join("results", record["SampleID"] + "." + bwa_prefix_name(bwa_prefix), record["Index"], record["FCID"])
    prefix = os.path.join(outdir, "reads")

    trim = PE_AdapterRemovalNode(input_files_1 = record["files"]["R1"], 
                                 input_files_2 = record["files"]["R2"], 
                                 output_prefix = os.path.join(outdir, "reads"))

    nodes = []
    for filename in ("aln", "unaln"):
        filename_prefix = prefix + ".singleton.%s.truncated" % (filename,)

        aln = SE_BWANode(input_file    = filename_prefix + ".gz",
                          output_file  = filename_prefix + ".bam",
                          prefix       = bwa_prefix,
                          read_group   = read_group,
                          dependencies = trim)
        val = ValidateBAMNode(config   = config, 
                              bamfile  = filename_prefix + ".bam",
                              ignore   = ["MATE_NOT_FOUND"],
                              dependencies = aln)
        nodes.append(val)


    aln  = PE_BWANode(input_file_1 = prefix + ".pair1.truncated.gz",
                        input_file_2 = prefix + ".pair2.truncated.gz",
                        output_file  = prefix + ".pairs.truncated.bam",
                        prefix       = bwa_prefix,
                        read_group   = read_group,
                        dependencies = trim)
    val = ValidateBAMNode(config   = config, 
                          bamfile  = prefix + ".pairs.truncated.bam",
                          ignore   = ["MATE_NOT_FOUND"],
                          dependencies = aln)
    nodes.append(val)


    return { "aligned"   : [prefix + ".singleton.aln.truncated.bam"],
             "unaligned" : [prefix + ".singleton.unaln.truncated.bam",
                            prefix + ".pairs.truncated.bam"],
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
        outdir = os.path.join("results", record["SampleID"] + "." + bwa_prefix_name(bwa_prefix))
        target = os.path.join(outdir, "%s.%s.bam" % (record["Index"], key))
        
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
            node = build_SE_nodes(bwa_prefix, record)
        else:
            assert False

        nodes[record["Index"]].append(node)
    
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
    parser.add_option("--run", action = "store_true", default = False)
    config, roots = parser.parse_args(argv)
    pipeline = pypeline.Pypeline(config)

    if not config.picard_root:
        ui.print_err("--picard-root must be set to the location of the Picard JARs.")
        return 1

    records = collect_records(roots)
    for bwa_prefix in config.bwa_prefix:
        for (sample, runs) in records.iteritems():
            pipeline.add_nodes(build_nodes(config, bwa_prefix, sample, runs))

    pipeline.run(dry_run = not config.run, max_running = 1)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

