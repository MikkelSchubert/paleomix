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
import sys

from pypeline import Pypeline
from pypeline.nodes.formats import FastaToPartitionedInterleavedPhyNode
from pypeline.nodes.raxml import *
from pypeline.nodes.mafft import MetaMAFFTNode

import common


def build_supermatrix(options, settings, destination, intervals, taxa, filtering_postfix, dependencies):
    input_files = {}
    for interval in intervals.itervalues():
        sequencedir = os.path.join(options.destination, "alignments", interval["Name"] + filtering_postfix)

        for sequence in common.collect_sequences(options, interval, taxa):
            filename = os.path.join(sequencedir, sequence + ".afa")
            record = {"name" : sequence}
            if interval["Protein coding"]:
                record["partition_by"] = ("12", "12", "3")

            assert filename not in input_files, filename
            input_files[filename] = record


    matrixprefix = os.path.join(destination, "alignments")
    supermatrix  = FastaToPartitionedInterleavedPhyNode(infiles        = input_files,
                                                        out_prefix     = matrixprefix,
                                                        partition_by   = "111",
                                                        dependencies   = dependencies)

    return RAxMLReduceNode(input_alignment  = matrixprefix + ".phy",
                           output_alignment = matrixprefix + ".reduced.phy",
                           input_partition  = matrixprefix + ".partitions",
                           output_partition = matrixprefix + ".reduced.partitions",
                           dependencies     = supermatrix)


def build_examl_nodes(options, settings, intervals, taxa, filtering, dependencies):
    filtering_postfix = ".filtered" if any(filtering.itervalues()) else ""
    destination = os.path.join(options.destination, "phylogenies", "examl.supermatrix" + filtering_postfix)

    input_alignment = os.path.join(destination, "alignments.reduced.phy")
    input_binary    = os.path.join(destination, "alignments.reduced.binary")

    supermatrix = build_supermatrix(options, settings, destination, intervals, taxa, filtering_postfix, dependencies)
    binary      = EXaMLParserNode(input_alignment = input_alignment,
                                  input_partition = os.path.join(destination, "alignments.reduced.partitions"),
                                  output_file     = input_binary,
                                  dependencies    = supermatrix)

    replicates = []
    for replicate_num in range(settings["ExaML"]["Replicates"]):
        replicate_destination = os.path.join(destination, "replicate_%04i" % replicate_num)
        initial_tree = os.path.join(replicate_destination, "initial.tree")

        tree = ParsimonatorNode(input_alignment = input_alignment,
                                output_tree     = initial_tree,
                                dependencies    = supermatrix)

        node = EXaMLNode(input_binary    = input_binary,
                         initial_tree    = initial_tree,
                         output_template = os.path.join(replicate_destination, "RAxML_%s"),
                         dependencies    = tree)

        replicates.append(node)


    return replicates


def chain_examl(pipeline, options, makefiles):
    destination = options.destination # Move to makefile
    for makefile in makefiles:
        settings  = makefile["Phylogenetic Inference"]
        taxa      = makefile["Project"]["Taxa"]
        intervals = makefile["Project"]["Intervals"]
        filtering = makefile["Project"]["Filter Singletons"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        makefile["Nodes"] = build_examl_nodes(options, settings, intervals, taxa, filtering, makefile["Nodes"])
    options.destination = destination
