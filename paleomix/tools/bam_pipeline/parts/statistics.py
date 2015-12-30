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

from paleomix.common.fileutils import \
    swap_ext

from paleomix.nodes.commands import \
    CoverageNode, \
    MergeCoverageNode, \
    DepthHistogramNode
from paleomix.tools.bam_pipeline.parts.summary import \
    SummaryTableNode


def add_statistics_nodes(config, makefile, target):
    features = makefile["Options"]["Features"]

    nodes = []
    if features["Depths"]:
        nodes.extend(_build_depth(config, target))

    if features["Summary"] or features["Coverage"]:
        make_summary = features["Summary"]
        coverage = _build_coverage(config, target, make_summary)
        if make_summary:
            summary_node = _build_summary_node(config, makefile,
                                               target, coverage)
            nodes.append(summary_node)
        elif features["Coverage"]:
            nodes.extend(coverage["Nodes"])

    target.nodes.extend(nodes)


def _build_summary_node(config, makefile, target, coverage):
    coverage_by_label = _build_coverage_nodes(config, target, use_label=True)

    return SummaryTableNode(config=config,
                            makefile=makefile,
                            target=target,
                            cov_for_lanes=coverage_by_label["Lanes"],
                            cov_for_libs=coverage_by_label["Libraries"],
                            dependencies=coverage["Nodes"])


def _build_depth(config, target):
    nodes = []
    for prefix in target.prefixes:
        for (roi_name, roi_filename) in _get_roi(prefix, name_prefix="."):
            if roi_filename is not None:
                # ROIs require indexed access, and hence that the final BAM
                # (either raw or realigned) has been built. By default, the
                # the realigned BAM is used (based on lexical order).
                bam_files = tuple(sorted(prefix.bams.items()))
                input_files, dependencies = bam_files[-1]
            else:
                input_files = {}
                for sample in prefix.samples:
                    input_files.update(sample.bams)
                dependencies = input_files.values()
                input_files = input_files.keys()

            output_filename = "%s.%s%s.depths" % (target.name, prefix.name,
                                                  roi_name)
            output_fpath = os.path.join(config.destination, output_filename)

            node = DepthHistogramNode(config=config,
                                      target_name=target.name,
                                      input_files=input_files,
                                      regions_file=roi_filename,
                                      output_file=output_fpath,
                                      dependencies=dependencies)
            nodes.append(node)

    return nodes


def _aggregate_for_prefix(cov, prefix, roi_name=None, into=None):
    prefix = _get_prefix_label(prefix, roi_name)
    results = {} if into is None else into
    for (key, files_and_nodes) in cov.iteritems():
        if prefix is None or (key[0] == prefix):
            results.update(files_and_nodes)
    return results


def _build_coverage(config, target, make_summary):
    merged_nodes = []
    coverage = _build_coverage_nodes(config, target)
    for prefix in target.prefixes:
        for (roi_name, _) in _get_roi(prefix):
            label = _get_prefix_label(prefix.name, roi_name)
            if not roi_name:
                postfix = prefix.name
            else:
                postfix = "%s.%s" % (prefix.name, roi_name)

            files_and_nodes = _aggregate_for_prefix(coverage["Libraries"],
                                                    label)
            output_filename = os.path.join(config.destination,
                                           "%s.%s.coverage"
                                           % (target.name, postfix))
            merged = MergeCoverageNode(input_files=files_and_nodes.keys(),
                                       output_file=output_filename,
                                       dependencies=files_and_nodes.values())

            merged_nodes.append(merged)

    files_and_nodes = _aggregate_for_prefix(coverage["Libraries"], None)
    if make_summary:
        files_and_nodes = _aggregate_for_prefix(coverage["Lanes"], None,
                                                into=files_and_nodes)

    all_nodes = []
    all_nodes.extend(files_and_nodes.itervalues())
    all_nodes.extend(merged_nodes)

    coverage["Nodes"] = tuple(all_nodes)

    return coverage


def _build_coverage_nodes(config, target, use_label=False):
    coverage = {"Lanes": collections.defaultdict(dict),
                "Libraries": collections.defaultdict(dict)}

    cache = {}
    for prefix in target.prefixes:
        for (roi_name, roi_filename) in _get_roi(prefix):
            prefix_label = prefix.label if use_label else prefix.name
            prefix_label = _get_prefix_label(prefix_label, roi_name)

            for sample in prefix.samples:
                for library in sample.libraries:
                    key = (prefix_label, target.name,
                           sample.name, library.name)

                    for lane in library.lanes:
                        for bams in lane.bams.values():
                            bams = _build_coverage_nodes_cached(config, bams,
                                                                target.name,
                                                                roi_name,
                                                                roi_filename,
                                                                cache)

                            coverage["Lanes"][key].update(bams)

                    bams = _build_coverage_nodes_cached(config, library.bams,
                                                        target.name, roi_name,
                                                        roi_filename, cache)
                    coverage["Libraries"][key].update(bams)
    return coverage


def _build_coverage_nodes_cached(config, files_and_nodes, target_name,
                                 roi_name, roi_filename, cache):
    output_ext = ".coverage"
    if roi_name:
        output_ext = ".%s.coverage" % roi_name

    coverages = {}
    for (input_filename, node) in files_and_nodes.iteritems():
        output_filename = swap_ext(input_filename, output_ext)

        cache_key = (roi_filename, input_filename)
        if cache_key not in cache:
            cache[cache_key] = CoverageNode(config=config,
                                            input_file=input_filename,
                                            output_file=output_filename,
                                            target_name=target_name,
                                            regions_file=roi_filename,
                                            dependencies=node)

        coverages[output_filename] = cache[cache_key]
    return coverages


def _get_roi(prefix, name_prefix=""):
    roi = [("", None)]
    for (name, path) in prefix.roi.iteritems():
        roi.append((name_prefix + name, path))
    return roi


def _get_prefix_label(label, roi_name):
    if not roi_name:
        return label
    return "%s:%s" % (label, roi_name)
