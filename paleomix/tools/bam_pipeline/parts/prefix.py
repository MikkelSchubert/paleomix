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

from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.nodes.picard import MergeSamFilesNode
from paleomix.tools.bam_pipeline.nodes import \
    index_and_validate_bam
from paleomix.nodes.validation import \
    DetectInputDuplicationNode

import paleomix.nodes.gatk as gatk


class Prefix:
    def __init__(self, config, prefix, samples, features, target):
        self.name      = prefix["Name"]
        self.label     = prefix.get("Label") or self.name
        self.reference = prefix["Reference"]
        self.roi       = prefix.get("RegionsOfInterest", {})

        self.samples = safe_coerce_to_tuple(samples)
        self.folder  = config.destination
        self.target  = target

        files_and_nodes = {}
        for sample in self.samples:
            files_and_nodes.update(sample.bams.iteritems())

        self.datadup_check = self._build_dataduplication_node(prefix, files_and_nodes)

        self.bams = {}
        if features["RawBAM"]:
            self.bams.update(self._build_raw_bam(config, prefix, files_and_nodes))
        if features["RealignedBAM"]:
            self.bams.update(self._build_realigned_bam(config, prefix, files_and_nodes))

        if not self.bams:
            for sample in self.samples:
                self.bams.update(sample.bams)

        nodes = [self.datadup_check]
        for sample in self.samples:
            nodes.extend(sample.nodes)
        self.nodes = tuple(nodes)

    def _build_raw_bam(self, config, prefix, files_and_bams):
        output_filename = os.path.join(self.folder, "%s.%s.bam" % (self.target, prefix["Name"]))
        validated_filename = os.path.join(self.folder, self.target, prefix["Name"] + ".validated")

        node = MergeSamFilesNode(config       = config,
                                 input_bams   = files_and_bams.keys(),
                                 output_bam   = output_filename,
                                 dependencies = self.datadup_check)
        validated_node = index_and_validate_bam(config, prefix, node, validated_filename)

        return {output_filename : validated_node}

    def _build_realigned_bam(self, config, prefix, bams):
        output_filename = os.path.join(self.folder, "%s.%s.realigned.bam" % (self.target, prefix["Name"]))
        intervals_filename = os.path.join(self.folder, self.target, prefix["Name"] + ".intervals")
        validated_filename = os.path.join(self.folder, self.target, prefix["Name"] + ".realigned.validated")

        trainer = gatk.GATKIndelTrainerNode(config=config,
                                            reference=prefix["Reference"],
                                            infiles=bams.keys(),
                                            outfile=intervals_filename,
                                            threads=config.gatk_max_threads,
                                            dependencies=self.datadup_check)

        aligner = gatk.GATKIndelRealignerNode(config=config,
                                              reference=prefix["Reference"],
                                              infiles=bams.keys(),
                                              intervals=intervals_filename,
                                              outfile=output_filename,
                                              dependencies=trainer)

        validated_node = index_and_validate_bam(config, prefix, aligner, validated_filename)

        return {output_filename: validated_node}

    def _build_dataduplication_node(self, prefix, files_and_nodes):
        destination = os.path.join(self.folder, self.target, prefix["Name"] + ".duplications_checked")
        return DetectInputDuplicationNode(input_files = files_and_nodes.keys(),
                                          output_file = destination,
                                          dependencies = files_and_nodes.values())
