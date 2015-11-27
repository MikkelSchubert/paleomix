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
import random
import collections

from paleomix.nodes.formats import \
    FastaToPartitionedInterleavedPhyNode as ToPhylipNode
from paleomix.nodes.raxml import \
     RAxMLParsimonyTreeNode
from paleomix.nodes.phylip import \
     PHYLIPBootstrapNode
from paleomix.nodes.examl import \
     ExaMLNode, \
     ExaMLParserNode
from paleomix.nodes.newick import \
     NewickRerootNode, \
     NewickSupportNode
from paleomix.common.fileutils import \
     swap_ext, \
     add_postfix


def _build_supermatrix(destination, input_files, exclude_samples, subset_files,
                       dependencies):
    matrixprefix = os.path.join(destination, "alignments")

    # This supermatrix (and partitions) is not used by the pipeline, but is
    # built for reference purposes; may be made optional in the future.
    supermatrix = ToPhylipNode(infiles=input_files,
                               out_prefix=matrixprefix,
                               exclude_groups=exclude_samples,
                               dependencies=dependencies,
                               file_dependencies=subset_files)

    # Supermatrix with empty columns (all N, n, or -) stripped
    reduced = ToPhylipNode(reduce=True,
                           infiles=input_files,
                           out_prefix=matrixprefix + ".reduced",
                           exclude_groups=exclude_samples,
                           dependencies=dependencies,
                           file_dependencies=subset_files)

    return (supermatrix, reduced)


def _examl_nodes(options, settings, input_alignment, input_partitions, input_binary, output_template, dependencies):
    parsimony_tree = output_template % ("parsimony_tree",)
    tree = RAxMLParsimonyTreeNode(input_alignment  = input_alignment,
                                  input_partitions = input_partitions,
                                  output_tree      = parsimony_tree,
                                  dependencies     = dependencies)

    params = ExaMLNode.customize(input_binary     = input_binary,
                                 initial_tree     = parsimony_tree,
                                 output_template  = output_template,
                                 threads          = options.examl_max_threads,
                                 dependencies     = tree)

    params.command.set_option("-m", settings["ExaML"]["Model"].upper())

    return params.build_node()


def _build_rerooted_trees(nodes, reroot_on):
    filenames = []
    for node in nodes:
        for filename in node.output_files:
            if filename.endswith(".result"):
                filenames.append(filename)

    output_file = os.path.dirname(filenames[0]) + ".newick"

    return NewickRerootNode(tree_files=filenames,
                            output_file=output_file,
                            taxa=reroot_on,
                            dependencies=nodes)


def _build_examl_replicates(options, phylo, destination, input_alignment, input_partition, dependencies):
    input_binary = os.path.join(destination, "alignments.reduced.binary")
    binary       = ExaMLParserNode(input_alignment = input_alignment,
                                   input_partition = input_partition,
                                   output_file     = input_binary,
                                   dependencies    = dependencies)

    replicates = []
    for replicate_num in range(phylo["ExaML"]["Replicates"]):
        replicate_destination = os.path.join(destination, "replicates")
        replicate_template    = os.path.join(replicate_destination, "replicate.%04i.%%s" % (replicate_num,))
        replicates.append(_examl_nodes(options, phylo, input_alignment, input_partition, input_binary, replicate_template, binary))

    if replicates:
        return _build_rerooted_trees(replicates, phylo["RootTreesOn"])

    return None


def _build_examl_bootstraps(options, phylo, destination, input_alignment, input_partition, dependencies):
    bootstraps = []
    num_bootstraps = phylo["ExaML"]["Bootstraps"]
    bootstrap_destination = os.path.join(destination, "bootstraps")
    bootstrap_template    = os.path.join(bootstrap_destination, "bootstrap.%04i.phy")

    for bootstrap_num in xrange(num_bootstraps):
        bootstrap_alignment = bootstrap_template % (bootstrap_num,)
        bootstrap = PHYLIPBootstrapNode(input_alignment  = input_alignment,
                                        input_partition  = input_partition,
                                        output_alignment = bootstrap_alignment,
                                        seed             = random.randint(1, 2**32 - 1),
                                        dependencies     = dependencies)

        bootstrap_binary      = swap_ext(bootstrap_alignment, ".binary")
        bootstrap_final       = swap_ext(bootstrap_alignment, ".%s")
        bs_binary   = ExaMLParserNode(input_alignment = bootstrap_alignment,
                                      input_partition = input_partition,
                                      output_file     = bootstrap_binary,
                                      dependencies    = bootstrap)

        bootstraps.append(_examl_nodes(options          = options,
                                       settings         = phylo,
                                       input_alignment  = bootstrap_alignment,
                                       input_partitions = input_partition,
                                       input_binary     = bootstrap_binary,
                                       output_template  = bootstrap_final,
                                       dependencies     = bs_binary))

    if bootstraps:
        return _build_rerooted_trees(bootstraps, phylo["RootTreesOn"])

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
                             dependencies       = (bootstrap, replicate)),


def _build_examl_nodes(options, settings, destination, input_files, subset_files, dependencies):
    input_alignment = os.path.join(destination, "alignments.reduced.phy")
    input_partition = os.path.join(destination, "alignments.reduced.partitions")

    excluded    = settings["ExcludeSamples"]
    supermatrix = _build_supermatrix(destination, input_files, excluded, subset_files, dependencies)
    examl_args  = (options, settings, destination, input_alignment, input_partition, supermatrix)

    examl_replicates = _build_examl_replicates(*examl_args)
    examl_bootstraps = _build_examl_bootstraps(*examl_args)
    examl_dependencies = add_bootstrap_support(destination, examl_replicates, examl_bootstraps)

    return examl_dependencies


def _build_examl_per_gene_nodes(options, settings, run_dd, roi, destination, filtering, dependencies):
    regions            = settings["Project"]["Regions"][roi["Name"]]
    sequences          = regions["Sequences"][roi["SubsetRegions"]]
    subset_files       = regions["SubsetFiles"][roi["SubsetRegions"]]
    filtering_postfix  = ".filtered" if any(filtering.itervalues()) else ""
    sequence_dir       = os.path.join(options.destination, "alignments", roi["Name"] + filtering_postfix)
    msa_enabled        = settings["MultipleSequenceAlignment"][regions["Name"]]["Enabled"]
    fasta_extension    = ".afa" if msa_enabled else ".fasta"

    partitions         = roi["Partitions"] or "111"

    nodes              = []
    for sequence in sequences:
        seq_source      = os.path.join(sequence_dir, sequence + fasta_extension)
        seq_destination = os.path.join(destination, sequence)
        input_files     = {sequence : {"partitions" : partitions, "filenames" : [seq_source]}}
        nodes.extend(_build_examl_nodes(options, run_dd, seq_destination, input_files, subset_files, dependencies))

    return nodes


def _build_examl_regions_nodes(options, settings, run_dd, destination, filtering, dependencies):
    input_files  = collections.defaultdict(dict)
    subset_files = []
    for (roi_name, roi_dd) in run_dd["RegionsOfInterest"].iteritems():
        regions     = settings["Project"]["Regions"][roi_name]
        subset_key  = roi_dd.get("SubsetRegions")
        sequences   = regions["Sequences"][subset_key]
        subset_files.extend(regions["SubsetFiles"][subset_key])

        partitions  = roi_dd["Partitions"]
        filtering_postfix  = ".filtered" if any(filtering.itervalues()) else ""
        sequence_dir       = os.path.join(options.destination, "alignments", roi_name + filtering_postfix)
        msa_enabled        = settings["MultipleSequenceAlignment"][regions["Name"]]["Enabled"]
        fasta_extension    = ".afa" if msa_enabled else ".fasta"

        if partitions:
            for sequence in sequences:
                seq_source = os.path.join(sequence_dir, sequence + fasta_extension)

                input_files[sequence] = {"partitions": partitions,
                                         "filenames": [seq_source]}
        else:
            filenames = []
            for sequence in sequences:
                filenames.append(os.path.join(sequence_dir, sequence + fasta_extension))

            input_files[roi_name] = {"partitions": "1",
                                     "filenames": filenames}

    return _build_examl_nodes(options, run_dd, destination, dict(input_files), subset_files, dependencies)


def build_phylogeny_nodes(options, settings, filtering, dependencies):
    nodes = []
    for (run_name, run_dd) in settings["PhylogeneticInference"].iteritems():
        destination = os.path.join(options.destination, "phylogenies", run_name)

        if run_dd["PerGeneTrees"]:
            run_nodes = []
            for roi in run_dd["RegionsOfInterest"].itervalues():
                roi_destination = os.path.join(destination, roi["Name"])
                run_nodes.extend(_build_examl_per_gene_nodes(options, settings, run_dd, roi, roi_destination, filtering, dependencies))
            nodes.extend(run_nodes)
        else:
            nodes.extend(_build_examl_regions_nodes(options, settings, run_dd, destination, filtering, dependencies))

    return nodes


def chain_examl(_pipeline, options, makefiles):
    destination = options.destination  # Move to makefile
    for makefile in makefiles:
        filtering = makefile["Project"]["FilterSingletons"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        makefile["Nodes"] = build_phylogeny_nodes(options, makefile, filtering, makefile["Nodes"])
    options.destination = destination
