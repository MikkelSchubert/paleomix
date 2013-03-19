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
import collections

from pypeline.node import MetaNode
from pypeline.common.fileutils import swap_ext
from pypeline.nodes.coverage import CoverageNode, MergeCoverageNode
from summary import SummaryTableNode

def add_statistics_nodes(config, makefile, target):
    features = makefile["Options"]["Features"]
    make_summary = "Summary" in features

    if "Coverage" not in features and not make_summary:
        return

    coverage = _build_coverage(config, makefile, target, make_summary)
    if make_summary:
        summary_node = SummaryTableNode(config         = config,
                                        makefile       = makefile,
                                        target         = target,
                                        cov_for_lanes  = coverage["Lanes"],
                                        cov_for_libs   = coverage["Libraries"],
                                        dependencies   = coverage["Nodes"])
        target.add_extra_nodes("Statistics", summary_node)
    else:
        target.add_extra_nodes("Coverage", coverage["Nodes"])


def _aggregate_for_prefix(cov, prefix, into = None):
    results = {} if into is None else {}
    for (key, files_and_nodes) in cov.iteritems():
        if prefix is None or (key[0] == prefix):
            results.update(files_and_nodes)
    return results


def _build_coverage(config, makefile, target, make_summary):
    merged_nodes = []
    coverage = _build_coverage_nodes(target)
    for prefix in target.prefixes:
        files_and_nodes = _aggregate_for_prefix(coverage["Libraries"], prefix.label)
        output_filename = os.path.join(config.destination, "%s.%s.coverage" % (target.name, prefix.name))
        merged = MergeCoverageNode(input_files  = files_and_nodes.keys(),
                                   output_file  = output_filename,
                                   dependencies = files_and_nodes.values())
        merged_nodes.append(merged)

    description = "Libraries"
    files_and_nodes = _aggregate_for_prefix(coverage["Libraries"], None)
    if make_summary:
        description = "Lanes and libraries"
        files_and_nodes = _aggregate_for_prefix(coverage["Lanes"], None, files_and_nodes)

    coverage["Nodes"] = [MetaNode(description = description,
                                  subnodes    = files_and_nodes.values()),
                         MetaNode(description = "Final coverage",
                                  subnodes    = merged_nodes)]

    return coverage


def _build_coverage_nodes(target):
    coverage = {"Lanes"     : collections.defaultdict(dict),
                "Libraries" : collections.defaultdict(dict)}
    for prefix in target.prefixes:
        for sample in prefix.samples:
            for library in sample.libraries:
                key = (prefix.label, target.name, sample.name, library.name)
                for lane in library.lanes:
                    for bams in lane.bams.values():
                        coverage["Lanes"][key].update(bams)
                coverage["Libraries"][key].update(library.bams)

    cov_by_filename = {}
    for files_and_nodes in coverage.itervalues():
        _build_coverage_nodes_cached(target.name, files_and_nodes,  cache = cov_by_filename)

    return coverage


def _build_coverage_nodes_cached(target_name, files_and_nodes, cache):
    for (key, filenames) in files_and_nodes.items():
        coverages = {}
        for (input_filename, node) in filenames.iteritems():
            output_filename = swap_ext(input_filename, ".coverage")
            if input_filename in cache:
                coverage = cache[input_filename]
            else:
                coverage = CoverageNode(input_file   = input_filename,
                                        output_file  = output_filename,
                                        target_name  = target_name,
                                        dependencies = node)
            coverages[output_filename] = coverage
        files_and_nodes[key] = coverages

