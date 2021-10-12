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

from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.nodes.picard import MergeSamFilesNode
from paleomix.nodes.validation import DetectInputDuplicationNode
from paleomix.pipelines.bam.nodes import index_and_validate_bam


class Prefix:
    def __init__(self, config, prefix, samples, features, target):
        self.name = prefix["Name"]
        self.roi = prefix.get("RegionsOfInterest", {})

        self.samples = safe_coerce_to_tuple(samples)
        self.folder = config.destination
        self.target = target

        files_and_nodes = {}
        for sample in self.samples:
            files_and_nodes.update(sample.bams.items())

        self.datadup_check = self._build_dataduplication_node(prefix, files_and_nodes)
        self.bams = self._build_bam(config, prefix, files_and_nodes)

        nodes = [self.datadup_check]
        for sample in self.samples:
            nodes.extend(sample.nodes)
        self.nodes = tuple(nodes)

    def _build_bam(self, config, prefix, files_and_bams):
        output_filename = os.path.join(
            self.folder, "%s.%s.bam" % (self.target, prefix["Name"])
        )
        validated_filename = os.path.join(
            self.folder, self.target, prefix["Name"] + ".validated"
        )

        node = MergeSamFilesNode(
            config=config,
            input_bams=list(files_and_bams),
            output_bam=output_filename,
            dependencies=self.datadup_check,
        )
        validated_node = index_and_validate_bam(
            config=config,
            prefix=prefix,
            node=node,
            log_file=validated_filename,
            validation_levels=("full", "partial"),
        )

        return {output_filename: validated_node}

    def _build_dataduplication_node(self, prefix, files_and_nodes):
        filename = prefix["Name"] + ".duplications_checked"
        destination = os.path.join(self.folder, self.target, filename)
        dependencies = list(files_and_nodes.values())

        return DetectInputDuplicationNode(
            input_files=list(files_and_nodes),
            output_file=destination,
            dependencies=dependencies,
        )
