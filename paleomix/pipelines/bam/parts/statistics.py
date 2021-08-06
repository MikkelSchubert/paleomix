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
import logging
import os

from paleomix.nodes.commands import CoverageNode, DepthHistogramNode


def add_statistics_nodes(config, makefile, target):
    features = makefile["Options"]["Features"]

    if features["Depths"]:
        target.nodes.extend(_build_depth_nodes(config, target))

    if features["Coverage"]:
        target.nodes.extend(_build_coverage_nodes(config, target))

    if features["Summary"]:
        log = logging.getLogger(__name__)
        log.warning("Summary report generation is currently disabled pending rework.")


def _build_depth_nodes(config, target):
    for prefix in target.prefixes:
        ((input_file, dependencies),) = prefix.bams.items()

        output_filename = "%s.%s.depths" % (target.name, prefix.name)
        output_fpath = os.path.join(config.destination, output_filename)

        yield DepthHistogramNode(
            target_name=target.name,
            input_file=input_file,
            output_file=output_fpath,
            dependencies=dependencies,
        )


def _build_coverage_nodes(config, target):
    for prefix in target.prefixes:
        ((input_file, dependencies),) = prefix.bams.items()

        output_filename = "%s.%s.coverage" % (target.name, prefix.name)
        output_fpath = os.path.join(config.destination, output_filename)

        yield CoverageNode(
            target_name=target.name,
            input_file=input_file,
            output_file=output_fpath,
            dependencies=dependencies,
        )
