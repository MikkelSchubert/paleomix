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
import copy

from pypeline.common.fileutils import missing_files
from pypeline.common.makefile import MakefileError
from pypeline.nodes.adapterremoval import SE_AdapterRemovalNode, \
                                          PE_AdapterRemovalNode

import pypeline.tools.bam_pipeline.paths as paths


class Reads:
    def __init__(self, config, record, quality_offset):
        self.quality_offset = quality_offset
        self.files = {}
        self.stats = None
        self.nodes = ()

        tags = record["Tags"]
        self.folder = os.path.join(config.destination, tags["Target"], "reads",
                                   tags["SM"], tags["LB"], tags["PU_cur"])

        lane_type = record.get("Type")
        if lane_type == "Raw":
            self._init_raw_reads(config, record)
        elif lane_type == "Trimmed":
            self.files.update(record["Data"])
        else:
            assert False, "Unexpected data type in Reads(): %s" % (repr(lane_type))

        for name in record["Options"]["ExcludeReads"]:
            self.files.pop(name, None)

        if config.allow_missing_input_files and self.nodes:
            input_missing  = missing_files(self.nodes[0].input_files)
            output_missing = missing_files(self.nodes[0].output_files)
            if input_missing and not output_missing:
                self.nodes = ()


    def _init_raw_reads(self, config, record):
        output_prefix = os.path.join(self.folder, "reads")
        files = record["Data"]
        if ("SE" in files):
            command = SE_AdapterRemovalNode.customize(input_files   = files["SE"],
                                                      output_prefix = output_prefix)
            self.files["Single"] = output_prefix + ".truncated.gz"
        else:
            command = PE_AdapterRemovalNode.customize(input_files_1 = files["PE_1"],
                                                      input_files_2 = files["PE_2"],
                                                      output_prefix = output_prefix)
            self.files["Single"]    = output_prefix + ".singleton.unaln.truncated.gz"
            self.files["Paired"]    = output_prefix + ".pair{Pair}.truncated.gz"
            self.files["Collapsed"] = output_prefix + ".singleton.aln.truncated.gz"
        self.stats = output_prefix + ".settings"

        quality_offset = self.quality_offset # record["Options"]["QualityOffset"]
        if quality_offset == "Solexa":
            quality_offset = 64
        command.command.set_parameter("--qualitybase", quality_offset)

        self.nodes = (command.build_node(),)
