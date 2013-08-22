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
from pypeline.nodes.depthhist import DepthHistogramNode
from pypeline.tools.bam_pipeline.parts.summary import SummaryTableNode


def add_statistics_nodes(config, makefile, target):
    features = set(makefile["Options"]["Features"])
    if not features & set(("Coverage", "Depths", "Summary")):
        return

    nodes = []
    if "Depths" in features:
        nodes.append(_build_depth(config, target))

    if "Summary" in features or "Coverage" in features:
        coverage = _build_coverage(config, target, ("Summary" in features))
        if "Summary" in features:
            coverage_by_label = _build_coverage_nodes(target, use_label = True)
            summary_node = SummaryTableNode(config         = config,
                                            makefile       = makefile,
                                            target         = target,
                                            cov_for_lanes  = coverage_by_label["Lanes"],
                                            cov_for_libs   = coverage_by_label["Libraries"],
                                            dependencies   = coverage["Node"])
            nodes.append(summary_node)
        elif "Coverage" in features:
            nodes.append(coverage["Node"])

    target.add_extra_nodes("Statistics", nodes)


def _build_depth(config, target):
    nodes = []
    for prefix in target.prefixes:
        input_files = {}
        for sample in prefix.samples:
            input_files.update(sample.bams)

        for (aoi_name, aoi_filename) in _get_aoi(prefix, name_prefix = "."):
            output_filename = os.path.join(config.destination,
                                           "%s.%s%s.depths" % (target.name, prefix.name, aoi_name))

            node = DepthHistogramNode(config         = config,
                                      target_name    = target.name,
                                      input_files    = input_files.keys(),
                                      intervals_file = aoi_filename,
                                      output_file    = output_filename,
                                      dependencies   = input_files.values())
            nodes.append(node)

    return MetaNode(description = "DepthHistograms",
                    subnodes    = nodes)



def _aggregate_for_prefix(cov, prefix, aoi_name = None, into = None):
    prefix = _get_prefix_label(prefix, aoi_name)
    results = {} if into is None else into
    for (key, files_and_nodes) in cov.iteritems():
        if prefix is None or (key[0] == prefix):
            results.update(files_and_nodes)
    return results


def _build_coverage(config, target, make_summary):
    merged_nodes = []
    coverage = _build_coverage_nodes(target)
    for prefix in target.prefixes:
        for (aoi_name, _) in _get_aoi(prefix):
            label = _get_prefix_label(prefix.name, aoi_name)
            postfix = prefix.name if (not aoi_name) else ("%s.%s" % (prefix.name, aoi_name))

            files_and_nodes = _aggregate_for_prefix(coverage["Libraries"], label)
            output_filename = os.path.join(config.destination, "%s.%s.coverage" % (target.name, postfix))
            merged = MergeCoverageNode(input_files  = files_and_nodes.keys(),
                                       output_file  = output_filename,
                                       dependencies = files_and_nodes.values())

            merged_nodes.append(merged)

    description = "Libraries"
    files_and_nodes = _aggregate_for_prefix(coverage["Libraries"], None)
    if make_summary:
        description = "Lanes and libraries"
        files_and_nodes = _aggregate_for_prefix(coverage["Lanes"], None, into = files_and_nodes)

    partial_nodes = MetaNode(description = description,
                             subnodes    = files_and_nodes.values())
    final_nodes   = MetaNode(description = "Final coverage",
                             subnodes    = merged_nodes)

    coverage["Node"] = MetaNode(description  = "Coverage",
                                dependencies = (partial_nodes, final_nodes))

    return coverage


def _build_coverage_nodes(target, use_label = False):
    coverage = {"Lanes"     : collections.defaultdict(dict),
                "Libraries" : collections.defaultdict(dict)}

    cache = {}
    for prefix in target.prefixes:
        for (aoi_name, aoi_filename) in _get_aoi(prefix):
            prefix_label = prefix.label if use_label else prefix.name
            prefix_label = _get_prefix_label(prefix_label, aoi_name)

            for sample in prefix.samples:
                for library in sample.libraries:
                    key = (prefix_label, target.name, sample.name, library.name)

                    for lane in library.lanes:
                        for bams in lane.bams.values():
                            bams = _build_coverage_nodes_cached(bams, target.name, aoi_name, aoi_filename, cache)
                            coverage["Lanes"][key].update(bams)

                    bams = _build_coverage_nodes_cached(library.bams, target.name, aoi_name, aoi_filename, cache)
                    coverage["Libraries"][key].update(bams)
    return coverage


def _build_coverage_nodes_cached(files_and_nodes, target_name, aoi_name, aoi_filename, cache):
    output_ext = ".coverage"
    if aoi_name:
        output_ext = ".%s.coverage" % aoi_name

    coverages = {}
    for (input_filename, node) in files_and_nodes.iteritems():
        output_filename = swap_ext(input_filename, output_ext)

        cache_key = (aoi_filename, input_filename)
        if cache_key not in cache:
            cache[cache_key] = CoverageNode(input_file     = input_filename,
                                            output_file    = output_filename,
                                            target_name    = target_name,
                                            intervals_file = aoi_filename,
                                            dependencies   = node)

        coverages[output_filename] = cache[cache_key]
    return coverages


def _get_aoi(prefix, name_prefix = ""):
    aoi = [("", None)]
    for (name, path) in prefix.aoi.iteritems():
        aoi.append((name_prefix + name, path))
    return aoi


def _get_prefix_label(label, aoi_name):
    if not aoi_name:
        return label
    return "%s:%s" % (label, aoi_name)
