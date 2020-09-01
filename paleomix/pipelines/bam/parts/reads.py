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

from paleomix.nodes.adapterremoval import SE_AdapterRemovalNode, PE_AdapterRemovalNode
from paleomix.nodes.validation import ValidateFASTQFilesNode


class Reads:
    def __init__(self, config, record, quality_offset):
        self.quality_offset = quality_offset
        self.files = {}
        self.stats = None
        self.validation = ()
        self.nodes = ()

        tags = record["Tags"]
        self.folder = os.path.join(
            config.destination,
            tags["Target"],
            "reads",
            tags["SM"],
            tags["LB"],
            tags["PU"],
            tags["Folder"],
        )

        lane_type = record.get("Type")
        if lane_type == "Raw":
            self._init_raw_reads(config, record)
        elif lane_type == "Trimmed":
            self._init_pretrimmed_reads(record)
        else:
            assert False, "Unexpected data type in Reads(): %s" % (repr(lane_type))

        for name, value in record["Options"]["ExcludeReads"].items():
            if value:
                self.files.pop(name, None)

    def _init_pretrimmed_reads(self, record):
        self.files.update(record["Data"])

        nodes = []
        validation = []
        for key, templ in self.files.items():
            if key == "Paired":
                input_files = [
                    templ.format(Pair=1),
                    templ.format(Pair=2),
                ]
            else:
                input_files = [templ]

            output_file = os.path.join(self.folder, key.lower() + ".statistics")
            node = ValidateFASTQFilesNode(
                input_files=input_files,
                output_file=output_file,
                offset=self.quality_offset,
                collapsed=("Collapsed" in key),
            )

            nodes.append(node)
            validation.append(output_file)

        self.nodes = tuple(nodes)
        self.validation = tuple(validation)

    def _init_raw_reads(self, config, record):
        ar_options = dict(record["Options"]["AdapterRemoval"])
        # Setup of "--collapsed" is handled by the node itself
        collapse_reads = ar_options.pop("--collapse")
        collapse_reads = collapse_reads or collapse_reads is None

        output_quality = self.quality_offset
        if output_quality == "Solexa":
            output_quality = "64"

        ar_options["--qualitybase"] = self.quality_offset
        ar_options["--qualitybase-output"] = output_quality

        init_args = {
            "output_prefix": os.path.join(self.folder, "reads"),
            "threads": config.adapterremoval_max_threads,
            "options": ar_options,
        }

        output_tmpl = "{output_prefix}.%s.gz".format(**init_args)

        if "SE" in record["Data"]:
            self.files["Single"] = output_tmpl % ("truncated",)
            init_args["input_file"] = record["Data"]["SE"]
            command = SE_AdapterRemovalNode(**init_args)
        else:
            self.files["Singleton"] = output_tmpl % ("singleton.truncated",)
            self.files["Paired"] = output_tmpl % ("pair{Pair}.truncated",)

            if collapse_reads:
                self.files["Collapsed"] = output_tmpl % ("collapsed",)
                self.files["CollapsedTruncated"] = output_tmpl % (
                    "collapsed.truncated",
                )

            init_args["collapse"] = collapse_reads
            init_args["input_file_1"] = record["Data"]["PE_1"]
            init_args["input_file_2"] = record["Data"]["PE_2"]
            command = PE_AdapterRemovalNode(**init_args)

        self.stats = os.path.join(self.folder, "reads.settings")
        self.nodes = (command,)
