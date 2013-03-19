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

from pypeline.common.utilities import safe_coerce_to_tuple
from pypeline.node import MetaNode
from pypeline.nodes.picard import MarkDuplicatesNode
from pypeline.tools.bam_pipeline.nodes import FilterUniqueBAMNode, \
                                              IndexAndValidateBAMNode, \
                                              MapDamageRescaleNode


class Library:
    def __init__(self, config, prefix, lanes, name):
        self.name    = name
        self.lanes   = safe_coerce_to_tuple(lanes)
        self.options = lanes[0].options
        self.folder  = os.path.dirname(self.lanes[0].folder)
        self.is_rmdupped = self.options["PCRDuplicates"]
        self.is_rescaled = self.options["RescaleQualities"]
        assert all((self.folder  == os.path.dirname(lane.folder)) for lane in self.lanes)

        bams = self._collect_bams(self.lanes)
        if self.is_rmdupped:
            bams = self._remove_pcr_duplicates(config, prefix, bams)

        if self.is_rescaled:
            bams = self._rescale_quality_scores(config, prefix, bams)

        self.bams = {}
        for files_and_nodes in bams.itervalues():
            self.bams.update(files_and_nodes)

        self.node = MetaNode(description  = "Library: %s" % os.path.basename(self.folder),
                             dependencies = self.bams.values())


    def _collect_bams(self, lanes):
        bams = {}
        for lane in lanes:
            for key, files in lane.bams.iteritems():
                key = "collapsed" if (key == "Collapsed") else "normal"
                bams.setdefault(key, {}).update(files)
        return bams


    def _remove_pcr_duplicates(self, config, prefix, bams):
        rmdup_cls = {"collapsed"  : FilterUniqueBAMNode,
                     "normal"     : MarkDuplicatesNode}

        results = {}
        for (key, files_and_nodes) in bams.items():
            output_filename = self.folder + ".rmdup.%s.bam" % key
            node = rmdup_cls[key](config       = config,
                                  input_bams   = files_and_nodes.keys(),
                                  output_bam   = output_filename,
                                  dependencies = files_and_nodes.values())
            validated_node = IndexAndValidateBAMNode(config, prefix, node)

            results[key] = {output_filename : validated_node}
        return results


    def _rescale_quality_scores(self, config, prefix, bams):
        files_and_nodes = {}
        for dd in bams.itervalues():
            files_and_nodes.update(dd)
        if not files_and_nodes:
            return bams

        output_filename = self.folder + ".rescaled.bam"
        node = MapDamageRescaleNode(config       = config,
                                    reference    = prefix["Reference"],
                                    input_files  = files_and_nodes.keys(),
                                    output_file  = output_filename,
                                    dependencies = files_and_nodes.values())
        validated_node = IndexAndValidateBAMNode(config, prefix, node)

        return {"Rescaled" : {output_filename : validated_node}}
