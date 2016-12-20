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
import copy

import paleomix.tools.bam_pipeline.paths as paths

from paleomix.atomiccmd.builder import \
    apply_options

from paleomix.common.makefile import \
    MakefileError

from paleomix.nodes.bwa import \
    BWAAlgorithmNode, \
    BWABacktrack, \
    BWASampe, \
    BWASamse
from paleomix.nodes.bowtie2 import \
    Bowtie2Node

from paleomix.tools.bam_pipeline.parts import \
    Reads
from paleomix.tools.bam_pipeline.nodes import \
    CleanupBAMNode, \
    index_and_validate_bam

from paleomix.common.fileutils import \
    swap_ext


#
_TRIMMED_READS_CACHE = {}


class Lane:
    def __init__(self, config, prefix, record, name):
        self.name = name
        self.bams = {}
        self.reads = None
        self.options = copy.deepcopy(record["Options"])
        self.tags = tags = copy.deepcopy(record["Tags"])
        self.tags["PU"] = self.tags["PU_src"]
        self.tags["PG"] = self.tags["PG"].lower()
        self.folder = os.path.join(config.destination,
                                   tags["Target"],
                                   prefix["Name"],
                                   tags["SM"],
                                   tags["LB"],
                                   tags["PU_cur"])

        if record["Type"] == "BAMs":
            self._init_pre_aligned_lane(config, prefix, record)
        else:
            self._init_reads(config, record)
            self._init_unaligned_lane(config, prefix, record)

    def _init_pre_aligned_lane(self, config, prefix, record):
        if prefix["Name"] not in record["Data"]:
            return

        input_filename = record["Data"][prefix["Name"]]
        output_filename = os.path.join(self.folder, "processed.bam")

        node = CleanupBAMNode(config=config,
                              reference=prefix["Reference"],
                              input_bam=input_filename,
                              output_bam=output_filename,
                              tags=self.tags,
                              dependencies=prefix["Nodes"])

        index_required = self._is_indexing_required(prefix)
        validated_node = index_and_validate_bam(config, prefix, node,
                                                create_index=index_required)

        self.bams["Processed"] = {output_filename: validated_node}

    def _init_reads(self, config, record):
        key = tuple(self.tags[key] for key in ("Target", "SM", "LB", "PU_cur"))
        if key not in _TRIMMED_READS_CACHE:
            _TRIMMED_READS_CACHE[key] \
                = Reads(config, record, record["Options"]["QualityOffset"])
        self.reads = _TRIMMED_READS_CACHE[key]

    def _init_unaligned_lane(self, config, prefix, record):
        prefix_key = "Nodes:%s" % (self.options["Aligners"]["Program"],)

        for (key, input_filename) in self.reads.files.iteritems():
            # Common parameters between BWA / Bowtie2
            output_filename = os.path.join(self.folder,
                                           "%s.bam" % (key.lower(),))

            parameters = {
                "input_file": input_filename,
                "output_file": output_filename,
                "prefix": prefix["Path"],
                "reference": prefix["Reference"],
                "dependencies": self.reads.nodes + prefix[prefix_key],
            }

            nodes = self._build_alignment_nodes(config=config,
                                                record=record,
                                                prefix=prefix,
                                                parameters=parameters)

            self.bams[key] = {output_filename: nodes}

    def _build_alignment_nodes(self, config, record, prefix, parameters):
        if self.options["Aligners"]["Program"] == "BWA":
            algorithm = self.options["Aligners"]["BWA"]["Algorithm"].lower()
            parameters["threads"] = config.bwa_max_threads
            parameters["algorithm"] = algorithm

            if algorithm == "backtrack":
                return self._build_bwa_backtrack(config=config,
                                                 record=record,
                                                 prefix=prefix,
                                                 parameters=parameters)
            elif algorithm in ('mem', 'bwasw'):
                return self._build_bwa_algorithm(config=config,
                                                 record=record,
                                                 prefix=prefix,
                                                 parameters=parameters)
            else:
                raise NotImplementedError('BWA %r not implemented!'
                                          % (algorithm,))

        elif self.options["Aligners"]["Program"] == "Bowtie2":
            return self._build_bowtie2(config=config,
                                       record=record,
                                       prefix=prefix,
                                       parameters=parameters)
        else:
            raise NotImplementedError('Aligner %r not implemented!'
                                      % (self.options["Aligners"]["Program"],))

    def _is_indexing_required(self, prefix):
        """Returns true if indexing lane BAMs is nessesary.
        """
        # Indexes are required for all files when calculating region statistics
        return bool(prefix.get("RegionsOfInterest")) or \
            (self.options["Features"]["RealignedBAM"] and not
             # and for BAMs fed to GATK, but in this case we only use these
             # indexes if we don't generate PCR filtered or recaled BAMs
             (self.options["Features"]["mapDamage"] == 'rescale' or
              self.options["PCRDuplicates"]))

    def _set_rg_tags(self, command):
        command.set_option("--rg-id", self.tags["ID"])
        for tag_name in ("SM", "LB", "PU", "PL", "PG"):
            tag_value = "%s:%s" % (tag_name, self.tags[tag_name])
            command.add_option("--rg", tag_value)

    ###########################################################################
    ###########################################################################
    # Construction of mapping nodes

    def _build_bwa_backtrack(self, config, prefix, record, parameters):
        if paths.is_paired_end(parameters["input_file"]):
            return self._build_bwa_backtrack_pe(config=config,
                                                prefix=prefix,
                                                record=record,
                                                parameters=parameters)
        else:
            return self._build_bwa_backtrack_se(config=config,
                                                prefix=prefix,
                                                record=record,
                                                parameters=parameters)

    def _build_bwa_backtrack_aln(self, parameters, input_file, output_file):
        """
        """
        node = BWABacktrack.customize(input_file=input_file,
                                      output_file=output_file,
                                      threads=parameters["threads"],
                                      prefix=parameters["prefix"],
                                      reference=parameters["reference"],
                                      dependencies=parameters["dependencies"])

        if not self.options["Aligners"]["BWA"]["UseSeed"]:
            node.commands["aln"].set_option("-l", 2 ** 16 - 1)

        if self.options["QualityOffset"] in (64, "Solexa"):
            node.commands["aln"].set_option("-I")

        apply_options(node.commands["aln"], self.options["Aligners"]["BWA"])

        return node.build_node()

    def _build_bwa_backtrack_se(self, config, prefix, record, parameters):
        input_file_fq = parameters.pop("input_file")
        output_file_bam = parameters.pop("output_file")
        output_file_sai = swap_ext(output_file_bam, ".sai")

        aln_node = self._build_bwa_backtrack_aln(parameters=parameters,
                                                 input_file=input_file_fq,
                                                 output_file=output_file_sai)

        sam_node = BWASamse.customize(input_file_fq=input_file_fq,
                                      input_file_sai=output_file_sai,
                                      output_file=output_file_bam,
                                      prefix=parameters["prefix"],
                                      reference=parameters["reference"],
                                      dependencies=aln_node)

        return self._finalize_nodes(config, prefix, parameters, sam_node)

    def _build_bwa_backtrack_pe(self, config, prefix, record, parameters):
        template = parameters.pop("input_file")
        output_bam = parameters.pop("output_file")

        aln_files = []
        aln_nodes = []
        for mate in (1, 2):
            input_file = template.format(Pair=mate)
            output_sai = swap_ext(output_bam, "%i.sai" % (mate,))

            aln_node = self._build_bwa_backtrack_aln(parameters=parameters,
                                                     input_file=input_file,
                                                     output_file=output_sai)

            aln_files.append(output_sai)
            aln_nodes.append(aln_node)

        sam_node = BWASampe.customize(input_file_sai_1=aln_files[0],
                                      input_file_sai_2=aln_files[1],
                                      input_file_fq_1=template.format(Pair=1),
                                      input_file_fq_2=template.format(Pair=2),
                                      output_file=output_bam,
                                      prefix=parameters['prefix'],
                                      reference=parameters["reference"],
                                      dependencies=aln_nodes)

        return self._finalize_nodes(config, prefix, parameters, sam_node)

    def _build_bwa_algorithm(self, config, prefix, record, parameters):
        if self.options["QualityOffset"] != 33:
            raise MakefileError("Mapping with BWA using the %r algorithm "
                                "currently does not support QualityOffsets "
                                "other than 33; please convert your FASTQ "
                                "if you wish to proceed.")

        self._set_pe_input_files(parameters)
        node = BWAAlgorithmNode.customize(**parameters)

        return self._finalize_nodes(config, prefix, parameters, node)

    def _build_bowtie2(self, config, prefix, record, parameters):
        self._set_pe_input_files(parameters)
        node = Bowtie2Node.customize(threads=config.bowtie2_max_threads,
                                     **parameters)

        command = node.commands["aln"]
        if self.options["QualityOffset"] == 33:
            command.set_option("--phred33")
        else:
            command.set_option("--phred64")

        for (key, value) in self.options["Aligners"]["Bowtie2"].iteritems():
            if key.startswith("-"):
                command.set_option(key, value)

        return self._finalize_nodes(config, prefix, parameters, node)

    def _finalize_nodes(self, config, prefix, parameters, node):
        self._set_rg_tags(node.commands["convert"])

        min_quality = self.options["Aligners"]["BWA"]["MinQuality"]
        node.commands["convert"].set_option('-q', min_quality)

        if self.options["Aligners"]["BWA"]["FilterUnmappedReads"]:
            node.commands["convert"].set_option('-F', "0x4")

        index_required = self._is_indexing_required(prefix)
        validated_node = index_and_validate_bam(config=config,
                                                prefix=parameters['prefix'],
                                                node=node.build_node(),
                                                create_index=index_required)

        return validated_node

    @classmethod
    def _set_pe_input_files(self, parameters):
        template = parameters.pop("input_file")
        if paths.is_paired_end(template):
            parameters["input_file_1"] = template.format(Pair=1)
            parameters["input_file_2"] = template.format(Pair=2)
        else:
            parameters["input_file_1"] = template
            parameters["input_file_2"] = None
