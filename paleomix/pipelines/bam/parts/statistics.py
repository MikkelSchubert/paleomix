#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
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

from paleomix.common.fileutils import swap_ext

from paleomix.nodes.commands import CoverageNode, MergeCoverageNode, DepthHistogramNode
from paleomix.pipelines.bam.parts.summary import SummaryTableNode


def add_statistics_nodes(config, makefile, target):
    features = makefile["Options"]["Features"]

    nodes = []
    if features["Depths"]:
        nodes.extend(_build_depth(config, target))

    if features["Summary"] or features["Coverage"]:
        make_summary = features["Summary"]
        coverage = _build_coverage(config, target, make_summary)
        if make_summary:
            summary_node = _build_summary_node(config, makefile, target, coverage)
            nodes.append(summary_node)
        elif features["Coverage"]:
            nodes.extend(coverage["Nodes"])

    target.nodes.extend(nodes)


def _build_summary_node(config, makefile, target, coverage):
    coverage_by_label = _build_coverage_nodes(target)

    return SummaryTableNode(
        config=config,
        makefile=makefile,
        target=target,
        cov_for_lanes=coverage_by_label["Lanes"],
        cov_for_libs=coverage_by_label["Libraries"],
        dependencies=coverage["Nodes"],
    )


def _build_depth(config, target):
    nodes = []
    for prefix in target.prefixes:
        ((input_file, dependencies),) = prefix.bams.items()

        output_filename = "%s.%s.depths" % (target.name, prefix.name)
        output_fpath = os.path.join(config.destination, output_filename)

        nodes.append(
            DepthHistogramNode(
                target_name=target.name,
                input_file=input_file,
                output_file=output_fpath,
                dependencies=dependencies,
            )
        )

    return nodes


def _aggregate_for_prefix(cov, prefix, into=None):
    results = {} if into is None else into
    for (key, files_and_nodes) in cov.items():
        if prefix is None or (key[0] == prefix):
            results.update(files_and_nodes)
    return results


def _build_coverage(config, target, make_summary):
    merged_nodes = []
    coverage = _build_coverage_nodes(target)
    for prefix in target.prefixes:
        files_and_nodes = _aggregate_for_prefix(coverage["Libraries"], prefix.name)
        output_filename = os.path.join(
            config.destination, "%s.%s.coverage" % (target.name, prefix.name)
        )
        merged = MergeCoverageNode(
            input_files=list(files_and_nodes.keys()),
            output_file=output_filename,
            dependencies=list(files_and_nodes.values()),
        )

        merged_nodes.append(merged)

    files_and_nodes = _aggregate_for_prefix(coverage["Libraries"], None)
    if make_summary:
        files_and_nodes = _aggregate_for_prefix(
            coverage["Lanes"], None, into=files_and_nodes
        )

    all_nodes = []
    all_nodes.extend(files_and_nodes.values())
    all_nodes.extend(merged_nodes)

    coverage["Nodes"] = tuple(all_nodes)

    return coverage


def _build_coverage_nodes(target):
    coverage = {
        "Lanes": collections.defaultdict(dict),
        "Libraries": collections.defaultdict(dict),
    }

    cache = {}
    for prefix in target.prefixes:
        for sample in prefix.samples:
            for library in sample.libraries:
                key = (prefix.name, target.name, sample.name, library.name)

                for lane in library.lanes:
                    for bams in lane.bams.values():
                        bams = _build_coverage_nodes_cached(bams, target.name, cache)

                        coverage["Lanes"][key].update(bams)

                bams = _build_coverage_nodes_cached(library.bams, target.name, cache)
                coverage["Libraries"][key].update(bams)
    return coverage


def _build_coverage_nodes_cached(files_and_nodes, target_name, cache):
    coverages = {}
    for (input_filename, node) in files_and_nodes.items():
        output_filename = swap_ext(input_filename, ".coverage")

        if input_filename not in cache:
            cache[input_filename] = CoverageNode(
                input_file=input_filename,
                output_file=output_filename,
                target_name=target_name,
                dependencies=node,
            )

        coverages[output_filename] = cache[input_filename]
    return coverages
