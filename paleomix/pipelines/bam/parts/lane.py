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
import copy

import paleomix.pipelines.bam.paths as paths

from paleomix.common.makefile import MakefileError

from paleomix.nodes.bwa import BWAAlgorithmNode, BWABacktrack, BWASampe, BWASamse
from paleomix.nodes.bowtie2 import Bowtie2Node

from paleomix.pipelines.bam.parts import Reads
from paleomix.pipelines.bam.nodes import index_and_validate_bam

from paleomix.common.fileutils import swap_ext


#
_TRIMMED_READS_CACHE = {}


class Lane:
    def __init__(self, config, prefix, record, name):
        self.name = name
        self.bams = {}
        self.reads = None
        self.options = copy.deepcopy(record["Options"])
        self.tags = tags = copy.deepcopy(record["Tags"])
        self.tags["PG"] = self.tags["PG"].lower()
        self.folder = os.path.join(
            config.destination,
            tags["Target"],
            prefix["Name"],
            tags["SM"],
            tags["LB"],
            tags["PU"],
            tags["Folder"],
        )

        self._init_reads(config, record)
        self._init_unaligned_lane(config, prefix, record)

    def _init_reads(self, config, record):
        key = tuple(self.tags[key] for key in ("Target", "SM", "LB", "PU", "DS"))
        if key not in _TRIMMED_READS_CACHE:
            _TRIMMED_READS_CACHE[key] = Reads(
                config, record, record["Options"]["QualityOffset"]
            )
        self.reads = _TRIMMED_READS_CACHE[key]

    def _init_unaligned_lane(self, config, prefix, record):
        prefix_key = "Nodes:%s" % (self.options["Aligners"]["Program"],)

        for (key, input_filename) in self.reads.files.items():
            # Common parameters between BWA / Bowtie2
            output_filename = os.path.join(self.folder, "%s.bam" % (key.lower(),))

            alignment_node = self._build_alignment_node(
                config=config,
                record=record,
                prefix=prefix,
                parameters={
                    "input_file": input_filename,
                    "output_file": output_filename,
                    "prefix": prefix["Path"],
                    "reference": prefix["Reference"],
                    "dependencies": self.reads.nodes + prefix[prefix_key],
                },
            )

            self.bams[key] = {
                output_filename: self._finalize_nodes(
                    config=config, prefix=prefix, node=alignment_node,
                )
            }

    def _build_alignment_node(self, config, record, prefix, parameters):
        if self.options["Aligners"]["Program"] == "BWA":
            algorithm = self.options["Aligners"]["BWA"]["Algorithm"].lower()
            parameters["threads"] = config.bwa_max_threads
            parameters["algorithm"] = algorithm

            if algorithm == "backtrack":
                return self._build_bwa_backtrack(
                    config=config, record=record, prefix=prefix, parameters=parameters
                )
            elif algorithm in ("mem", "bwasw"):
                return self._build_bwa_algorithm(
                    config=config, record=record, prefix=prefix, parameters=parameters
                )
            else:
                raise NotImplementedError("BWA %r not implemented!" % (algorithm,))

        elif self.options["Aligners"]["Program"] == "Bowtie2":
            return self._build_bowtie2(
                config=config, record=record, prefix=prefix, parameters=parameters
            )
        else:
            raise NotImplementedError(
                "Aligner %r not implemented!" % (self.options["Aligners"]["Program"],)
            )

    ###########################################################################
    ###########################################################################
    # Construction of mapping nodes

    def _build_bwa_backtrack(self, config, prefix, record, parameters):
        if paths.is_paired_end(parameters["input_file"]):
            return self._build_bwa_backtrack_pe(
                config=config, prefix=prefix, record=record, parameters=parameters
            )
        else:
            return self._build_bwa_backtrack_se(
                config=config, prefix=prefix, record=record, parameters=parameters
            )

    def _build_bwa_backtrack_aln(self, parameters, input_file, output_file):
        options = dict(self.options["Aligners"]["BWA"])
        if not self.options["Aligners"]["BWA"]["UseSeed"]:
            options["-l"] = 2 ** 16 - 1

        if self.options["QualityOffset"] in (64, "Solexa"):
            options["-I"] = None

        return BWABacktrack(
            input_file=input_file,
            output_file=output_file,
            threads=parameters["threads"],
            prefix=parameters["prefix"],
            reference=parameters["reference"],
            mapping_options=options,
            dependencies=parameters["dependencies"],
        )

    def _build_bwa_backtrack_se(self, config, prefix, record, parameters):
        input_file_fq = parameters.pop("input_file")
        output_file_bam = parameters.pop("output_file")
        output_file_sai = swap_ext(output_file_bam, ".sai")

        aln_node = self._build_bwa_backtrack_aln(
            parameters=parameters, input_file=input_file_fq, output_file=output_file_sai
        )

        return BWASamse(
            input_file_fq=input_file_fq,
            input_file_sai=output_file_sai,
            output_file=output_file_bam,
            prefix=parameters["prefix"],
            reference=parameters["reference"],
            mapping_options=self.options["Aligners"]["BWA"],
            cleanup_options=self._cleanup_options("BWA"),
            dependencies=aln_node,
        )

    def _build_bwa_backtrack_pe(self, config, prefix, record, parameters):
        template = parameters.pop("input_file")
        output_bam = parameters.pop("output_file")

        output_sai_1 = swap_ext(output_bam, "%i.sai" % (1,))
        aln_node_1 = self._build_bwa_backtrack_aln(
            parameters=parameters,
            input_file=template.format(Pair=1),
            output_file=output_sai_1,
        )

        output_sai_2 = swap_ext(output_bam, "%i.sai" % (2,))
        aln_node_2 = self._build_bwa_backtrack_aln(
            parameters=parameters,
            input_file=template.format(Pair=2),
            output_file=output_sai_2,
        )

        return BWASampe(
            input_file_sai_1=output_sai_1,
            input_file_sai_2=output_sai_2,
            input_file_fq_1=template.format(Pair=1),
            input_file_fq_2=template.format(Pair=2),
            output_file=output_bam,
            prefix=parameters["prefix"],
            reference=parameters["reference"],
            mapping_options=self.options["Aligners"]["BWA"],
            cleanup_options=self._cleanup_options("BWA"),
            dependencies=(aln_node_1, aln_node_2),
        )

    def _build_bwa_algorithm(self, config, prefix, record, parameters):
        if self.options["QualityOffset"] != 33:
            raise MakefileError(
                "Mapping with BWA using the %r algorithm currently does not support "
                "QualityOffsets other than 33; please convert your FASTQ if you wish "
                "to proceed."
            )

        parameters = self._set_pe_input_files(parameters)
        parameters["mapping_options"] = self.options["Aligners"]["BWA"]
        parameters["cleanup_options"] = self._cleanup_options("BWA")

        return BWAAlgorithmNode(**parameters)

    def _build_bowtie2(self, config, prefix, record, parameters):
        parameters = self._set_pe_input_files(parameters)
        parameters["threads"] = config.bowtie2_max_threads

        parameters["mapping_options"] = dict(self.options["Aligners"]["Bowtie2"])
        parameters["cleanup_options"] = self._cleanup_options("Bowtie2")

        if self.options["QualityOffset"] != 33:
            parameters["mapping_options"]["--phred64"] = None

        return Bowtie2Node(**parameters)

    def _cleanup_options(self, aligner):
        options = {
            "--rg-id": self.tags["ID"],
            "--rg": [
                "%s:%s" % (tag_name, self.tags[tag_name])
                for tag_name in ("SM", "LB", "PU", "PL", "PG", "DS")
            ],
            "-q": self.options["Aligners"][aligner]["MinQuality"],
        }

        if self.options["Aligners"][aligner]["FilterUnmappedReads"]:
            options["-F"] = "0x4"

        return options

    def _finalize_nodes(self, config, prefix, node):
        # Indexes are required for all files when calculating region statistics
        return index_and_validate_bam(
            config=config,
            prefix=prefix,
            node=node,
            create_index=bool(prefix.get("RegionsOfInterest")),
        )

    @staticmethod
    def _set_pe_input_files(parameters):
        parameters = dict(parameters)
        template = parameters.pop("input_file")
        if paths.is_paired_end(template):
            parameters["input_file_1"] = template.format(Pair=1)
            parameters["input_file_2"] = template.format(Pair=2)
        else:
            parameters["input_file_1"] = template
            parameters["input_file_2"] = None

        return parameters
