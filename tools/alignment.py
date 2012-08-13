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
from pypeline.common.utilities import safe_coerce_to_tuple


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


def build_full_path(record, bwa_prefix = ""):
    if bwa_prefix:
        bwa_prefix = "." + bwa_prefix_name(bwa_prefix)

    return os.path.join("results", record["Name"] + bwa_prefix, record["Sample"], record["Library"], record["Barcode"])


def build_sample_path(record, bwa_prefix = ""):
    if bwa_prefix:
        bwa_prefix = "." + bwa_prefix_name(bwa_prefix)

    return os.path.join("results", record["Name"] + bwa_prefix, record["Sample"])



def collect_files(record):
    """ """
    template = record["Path"]
    if template.format(Pair = 1) == template:
        # Single-ended only
        return {"R1" : list(sorted(glob.glob(template))),
                "R2" : []}
    
    files = {}
    for (ii, name) in enumerate(("R1", "R2"), start = 1):
        files[name] = list(sorted(glob.glob(template.format(Pair = ii))))
    return files


def validate_records_unique(records):
    paths = collections.defaultdict(list)
    for (name, runs) in records.iteritems():
        for record in runs:
            paths[build_full_path(record)].append(record)

    errors = False
    for (path, records) in paths.iteritems():
        if len(records) > 1:
            errors = True
            ui.print_err("ERROR: {0} too similar records found (combination of specified fields must be unique):".format(len(records)))
            ui.print_err("\t- Name:    {Name}\n\t\t- Sample:  {Sample}\n\t\t- Library: {Library}\n\t\t- Barcode: {Barcode}\n".format(**record))

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
    for (target, runs) in records.iteritems():
        for record in runs:
            template = record["Path"]
            current_paths = [template]
            if template.format(Pair = 1) != template:
                for end in (1, 2):
                    current_paths.append(template.format(Pair = end))
            
            for path in current_paths:
                paths[path].append(record)

    printed = list()
    for path in sorted(paths):
        current_records = tuple(sorted(paths[path]))
        if (len(current_records) > 1) and (current_records not in printed):
            printed.append(current_records)
            errors = True
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

    return not errors


def collect_records(filenames):
    def _split(line):
        return filter(None, line.strip().split("\t"))

    records = collections.defaultdict(list)
    for filename in filenames:
        with open(filename) as mkfile:
            header = _split(mkfile.readline())

            for line in mkfile:
                if line.startswith("#") or not line.strip():
                    continue

                record = dict(zip(header, _split(line)))
                record["files"] = collect_files(record)
                records[record["Name"]].append(record)
                
    return records


def validate_bams(config, node):
    validators = []
    for filename in node.output_files:
        if not filename.lower().endswith(".bam"):
            continue
                
        validators.append(ValidateBAMNode(config   = config, 
                                          bamfile  = filename,
                                          ignore   = ["MATE_NOT_FOUND"],
                                          dependencies = node))
    assert len(validators) == 1
    return validators[0]
            


_ADAPTERRM_SE_CACHE = {}
def build_SE_nodes(config, bwa_prefix, record):
    try:
        trim_prefix, trim_node = _ADAPTERRM_SE_CACHE[id(record)]
    except KeyError:
        trim_prefix = os.path.join(build_full_path(record), "reads")
        trim_node   = SE_AdapterRemovalNode(record["files"]["R1"], trim_prefix)
        _ADAPTERRM_SE_CACHE[id(record)] = (trim_prefix, trim_node)


    output_dir  = build_full_path(record, bwa_prefix)
    output_file = os.path.join(output_dir, "reads.truncated.bam")
    
    read_group = "@RG\tID:%(Library)s\tSM:%(Sample)s\tLB:%(Library)s\tPU:%(Barcode)s\tPL:Illumina" % record
    aln  = SE_BWANode(input_file   = trim_prefix + ".truncated.gz",
                      output_file  = output_file,
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      threads      = config.max_bwa_threads,
                      dependencies = trim_node)

    return { "aligned"   : {"files" : [], 
                            "nodes" : []},
             "unaligned" : {"files" : [output_file],
                            "nodes" : [validate_bams(config, aln)]},
             "record"    : record}


_ADAPTERRM_PE_CACHE = {}
def build_PE_nodes(config, bwa_prefix, record):
    try:
        trim_prefix, trim_node = _ADAPTERRM_PE_CACHE[id(record)]
    except KeyError:
        trim_prefix = os.path.join(build_full_path(record), "reads")
        trim_node   = PE_AdapterRemovalNode(input_files_1 = record["files"]["R1"], 
                                            input_files_2 = record["files"]["R2"], 
                                            output_prefix = trim_prefix)
        _ADAPTERRM_PE_CACHE[id(record)] = (trim_prefix, trim_node)


    nodes = {}
    output_dir = build_full_path(record, bwa_prefix)
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
        nodes[filename] = [validate_bams(config, aln)]


    output_filename = os.path.join(output_dir, "reads.pairs.truncated.bam")
    aln  = PE_BWANode(input_file_1 = trim_prefix + ".pair1.truncated.gz",
                      input_file_2 = trim_prefix + ".pair2.truncated.gz",
                      output_file  = output_filename,
                      prefix       = bwa_prefix,
                      read_group   = read_group,
                      threads      = config.max_bwa_threads,
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

        record = node["record"]
        outdir = build_sample_path(record, bwa_prefix)
        target = os.path.join(outdir, "%s.%s.bam" % (record["Library"], key))
        
        merge = MergeSamFilesNode(config       = config, 
                                  input_files  = files, 
                                  output_file  = target,
                                  create_index = False,
                                  dependencies = deps)

        yield (key, target, validate_bams(config, merge))

    
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

    return filename, validate_bams(config, dedup)

    


def build_nodes(config, bwa_prefix, sample, records):
    nodes = collections.defaultdict(list)
    for record in records:
        if record["files"]["R2"]:
            node = build_PE_nodes(config, bwa_prefix, record)
        elif record["files"]["R1"]:
            node = build_SE_nodes(config, bwa_prefix, record)
        else:
            ui.print_err("Could not find any files for record:\n\t- Name: %(Name)s\n\t- Sample: %(Sample)s\n\t- Library: %(Library)s\n\t- Barcode: %(Barcode)s" % record)
            return None

        nodes[record["Library"]].append(node)
    
    deduped = []
    for nodes in nodes.values():
        for (key, filename, node) in build_lib_merge_nodes(config, bwa_prefix, nodes):
            deduped.append(build_dedupe_node(config, key, filename, node))


    target = os.path.join("results", "%s.%s.bam" % (sample, bwa_prefix_name(bwa_prefix)))
    merge = MergeSamFilesNode(config       = config,
                              input_files  = [filename for (filename, _) in deduped],
                              output_file  = target,
                              dependencies = [node for (_, node) in deduped])

    return validate_bams(config, merge)
    


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
    if not validate_records(records):
        return 1

    for bwa_prefix in config.bwa_prefix:
        for (sample, runs) in records.iteritems():
            nodes = build_nodes(config, bwa_prefix, sample, runs)
            if not nodes:
                return 1

            pipeline.add_nodes(nodes)

    pipeline.run(dry_run = not config.run, max_running = config.max_threads)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

