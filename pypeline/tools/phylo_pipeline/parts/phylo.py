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
#!/usr/bin/python

import os

from pypeline.node import MetaNode
from pypeline.nodes.formats import \
     FastaToPartitionedInterleavedPhyNode
from pypeline.nodes.raxml import \
     RAxMLReduceNode, \
     RAxMLBootstrapNode
from pypeline.nodes.examl import \
     EXaMLNode, \
     EXaMLParserNode, \
     ParsimonatorNode
from pypeline.nodes.newick import \
     NewickRerootNode, \
     NewickSupportNode
from pypeline.common.fileutils import \
     swap_ext, \
     add_postfix

import pypeline.tools.phylo_pipeline.parts.common as common


# Number of files to run per bootstrap generation;
# chunking is used to reduce the overhead
_BOOTSTRAP_CHUNK = 25


def build_supermatrix(options, settings, afa_ext, destination, interval, filtering_postfix, dependencies):
    input_files = {}
    sequencedir = os.path.join(options.destination, "alignments", interval["Name"] + filtering_postfix)
    for sequence in common.get_sequences(options, interval):
        filename = os.path.join(sequencedir, sequence + afa_ext)
        record = {"name" : sequence}
        if interval["Protein coding"]:
            record["partition_by"] = ("12", "12", "3")

        assert filename not in input_files, filename
        input_files[filename] = record

    excluded_groups = settings["ExcludeGroups"]
    matrixprefix = os.path.join(destination, "alignments")
    supermatrix  = FastaToPartitionedInterleavedPhyNode(infiles        = input_files,
                                                        out_prefix     = matrixprefix,
                                                        partition_by   = "111",
                                                        exclude_groups = excluded_groups,
                                                        dependencies   = dependencies)

    return RAxMLReduceNode(input_alignment  = matrixprefix + ".phy",
                           output_alignment = matrixprefix + ".reduced.phy",
                           input_partition  = matrixprefix + ".partitions",
                           output_partition = matrixprefix + ".reduced.partitions",
                           dependencies     = supermatrix)


def _examl_nodes(options, settings, input_alignment, input_binary, output_template, dependencies):
    initial_tree = output_template % ("parsimony_tree",)
    tree = ParsimonatorNode(input_alignment = input_alignment,
                            output_tree     = initial_tree,
                            dependencies    = dependencies)

    params = EXaMLNode.customize(input_binary    = input_binary,
                                 initial_tree    = initial_tree,
                                 output_template = output_template,
                                 threads         = options.max_examl_threads,
                                 dependencies    = tree)

    params.command.set_option("-m", settings["ExaML"]["Model"].upper())

    return params.build_node()


def _build_rerooted_trees(meta_node):
    filenames = []
    for node in meta_node.subnodes:
        for filename in node.output_files:
            if filename.endswith(".result"):
                filenames.append(filename)

    output_file = os.path.dirname(filenames[0]) + ".newick"
    output_node = NewickRerootNode(tree_files   = filenames,
                                   output_file  = output_file,
                                   dependencies = meta_node)
    return output_node


def build_examl_replicates(options, phylo, destination, input_alignment, input_partition, dependencies):
    input_binary = os.path.join(destination, "alignments.reduced.binary")
    binary       = EXaMLParserNode(input_alignment = input_alignment,
                                   input_partition = input_partition,
                                   output_file     = input_binary,
                                   dependencies    = dependencies)

    replicates = []
    for replicate_num in range(phylo["ExaML"]["Replicates"]):
        replicate_destination = os.path.join(destination, "replicates")
        replicate_template    = os.path.join(replicate_destination, "replicate.%04i.%%s" % (replicate_num,))
        replicates.append(_examl_nodes(options, phylo, input_alignment, input_binary, replicate_template, binary))

    if replicates:
        meta = MetaNode(description  = "Replicates",
                        subnodes     = replicates,
                        dependencies = binary)

        return _build_rerooted_trees(meta)
    return None


def build_examl_bootstraps(options, phylo, destination, input_alignment, input_partition, dependencies):
    bootstraps = []
    num_bootstraps = phylo["ExaML"]["Bootstraps"]
    for bootstrap_start in range(0, num_bootstraps, _BOOTSTRAP_CHUNK):
        bootstrap_destination = os.path.join(destination, "bootstraps")
        bootstrap_template    = os.path.join(bootstrap_destination, "bootstrap.%04i.phy")

        bootstrap   = RAxMLBootstrapNode(input_alignment  = input_alignment,
                                         input_partition  = input_partition,
                                         template         = bootstrap_template,
                                         bootstraps       = _BOOTSTRAP_CHUNK,
                                         start            = bootstrap_start,
                                         dependencies     = dependencies)

        bootstrap_end = min(num_bootstraps, bootstrap_start + _BOOTSTRAP_CHUNK)
        for bootstrap_num in range(bootstrap_start, bootstrap_end):
            bootstrap_alignment   = bootstrap_template % (bootstrap_num,)
            bootstrap_binary      = swap_ext(bootstrap_alignment, ".binary")
            bootstrap_final       = swap_ext(bootstrap_alignment, ".%s")
            bs_binary   = EXaMLParserNode(input_alignment = bootstrap_alignment,
                                          input_partition = input_partition,
                                          output_file     = bootstrap_binary,
                                          dependencies    = bootstrap)

            bootstraps.append(_examl_nodes(options         = options,
                                           settings        = phylo,
                                           input_alignment = bootstrap_alignment,
                                           input_binary    = bootstrap_binary,
                                           output_template = bootstrap_final,
                                           dependencies    = bs_binary))

    if bootstraps:
        meta = MetaNode(description  = "Bootstraps",
                        subnodes     = bootstraps,
                        dependencies = dependencies)
        return _build_rerooted_trees(meta)
    return None


def add_bootstrap_support(destination, replicate, bootstrap):
    if not (replicate and bootstrap):
        return filter(None, (replicate, bootstrap))

    replicate_file = os.path.join(destination, "replicates.newick")
    bootstrap_file = os.path.join(destination, "bootstraps.newick")
    output_file    = add_postfix(replicate_file, ".support")

    return NewickSupportNode(main_tree_files    = replicate_file,
                             support_tree_files = bootstrap_file,
                             output_file        = output_file,
                             dependencies       = (bootstrap, replicate))


def build_examl_nodes(options, settings, intervals, filtering, dependencies):
    examl_nodes = []
    phylo = settings["Phylogenetic Inference"]
    filtering_postfix = ".filtered" if any(filtering.itervalues()) else ""
    for interval in intervals.itervalues():
        destination = os.path.join(options.destination, "phylogenies", "examl", interval["Name"] + filtering_postfix)
        afa_ext = ".afa" if settings["MSAlignment"]["Enabled"] else ".fasta"

        input_alignment = os.path.join(destination, "alignments.reduced.phy")
        input_partition = os.path.join(destination, "alignments.reduced.partitions")

        supermatrix = build_supermatrix(options, phylo, afa_ext, destination, interval, filtering_postfix, dependencies)
        examl_args  = (options, phylo, destination, input_alignment, input_partition, supermatrix)

        examl_replicates = build_examl_replicates(*examl_args)
        examl_bootstraps = build_examl_bootstraps(*examl_args)
        examl_dependencies = add_bootstrap_support(destination, examl_replicates, examl_bootstraps)

        if examl_dependencies:
            examl_nodes.append(MetaNode(description  = interval["Name"],
                                        dependencies = examl_dependencies))

    return MetaNode(description = "EXaML",
                    dependencies = examl_nodes)


def chain_examl(pipeline, options, makefiles):
    destination = options.destination # Move to makefile
    for makefile in makefiles:
        intervals = makefile["Project"]["Intervals"]
        filtering = makefile["Project"]["Filter Singletons"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        makefile["Nodes"] = build_examl_nodes(options, makefile, intervals, filtering, makefile["Nodes"])
    options.destination = destination
