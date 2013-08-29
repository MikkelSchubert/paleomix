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
     RAxMLBootstrapNode, \
     EXaMLNode, \
     EXaMLParserNode, \
     ParsimonatorNode
from pypeline.common.fileutils import swap_ext

import pypeline.tools.phylo_pipeline.parts.common as common


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


def _examl_nodes(settings, input_alignment, input_binary, output_template, dependencies):
    initial_tree = output_template % ("parsimony_tree",)
    tree = ParsimonatorNode(input_alignment = input_alignment,
                            output_tree     = initial_tree,
                            dependencies    = dependencies)

    return EXaMLNode(input_binary    = input_binary,
                     initial_tree    = initial_tree,
                     output_template = output_template,
                     threads         = settings["ExaML"]["Threads"],
                     dependencies    = tree)


def build_examl_replicates(phylo, destination, input_alignment, input_partition, dependencies):
    input_binary = os.path.join(destination, "alignments.reduced.binary")
    binary       = EXaMLParserNode(input_alignment = input_alignment,
                                   input_partition = input_partition,
                                   output_file     = input_binary,
                                   dependencies    = dependencies)

    replicates = []
    for replicate_num in range(phylo["ExaML"]["Replicates"]):
        replicate_destination = os.path.join(destination, "replicate_%04i" % replicate_num)
        replicate_template    = os.path.join(replicate_destination, "EXaML_%s")
        replicates.append(_examl_nodes(phylo, input_alignment, input_binary, replicate_template, binary))

    if replicates:
        return [MetaNode(description  = "Replicates",
                         subnodes     = replicates,
                         dependencies = binary)]

    return []


def build_examl_bootstraps(phylo, destination, input_alignment, input_partition, dependencies):
    bootstraps = []
    for bootstrap_start in range(0, phylo["ExaML"]["Bootstraps"], 50):
        bootstrap_destination = os.path.join(destination, "bootstraps")
        bootstrap_template    = os.path.join(bootstrap_destination, "bootstrap.%04i.phy")

        bootstrap   = RAxMLBootstrapNode(input_alignment  = input_alignment,
                                         input_partition  = input_partition,
                                         template         = bootstrap_template,
                                         bootstraps       = 50,
                                         start            = bootstrap_start,
                                         dependencies     = dependencies)

        for bootstrap_num in range(bootstrap_start, bootstrap_start + 50):
            bootstrap_alignment   = bootstrap_template % (bootstrap_num,)
            bootstrap_binary      = swap_ext(bootstrap_alignment, ".binary")
            bootstrap_final       = swap_ext(bootstrap_alignment, ".EXaML_%s")
            bs_binary   = EXaMLParserNode(input_alignment = bootstrap_alignment,
                                          input_partition = input_partition,
                                          output_file     = bootstrap_binary,
                                          dependencies    = bootstrap)

            bootstraps.append(_examl_nodes(settings        = phylo,
                                           input_alignment = bootstrap_alignment,
                                           input_binary    = bootstrap_binary,
                                           output_template = bootstrap_final,
                                           dependencies    = bs_binary))

    if bootstraps:
        return [MetaNode(description  = "Bootstraps",
                         subnodes     = bootstraps,
                         dependencies = dependencies)]
    return []


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
        examl_args  = (phylo, destination, input_alignment, input_partition, supermatrix)

        examl_dependencies = []
        examl_dependencies.extend(build_examl_replicates(*examl_args))
        examl_dependencies.extend(build_examl_bootstraps(*examl_args))

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
