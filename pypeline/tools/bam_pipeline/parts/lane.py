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

import reads
from pypeline.node import MetaNode
from pypeline.nodes.bwa import BWANode
from pypeline.nodes.bowtie2 import Bowtie2Node
from pypeline.common.utilities import safe_coerce_to_tuple

import pypeline.tools.bam_pipeline.paths as paths
from pypeline.tools.bam_pipeline.nodes import CleanupBAMNode, \
                                              IndexAndValidateBAMNode


#
_TRIMMED_READS_CACHE = {}


class Lane:
    def __init__(self, config, prefix, record, name):
        self.name    = name
        self.bams    = {}
        self.reads   = None
        self.options = copy.deepcopy(record["Options"])
        self.tags    = tags = copy.deepcopy(record["Tags"])
        self.folder  = os.path.join(config.destination, tags["Target"], prefix["Name"],
                                    tags["SM"], tags["LB"], tags["PU"])

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

        node = CleanupBAMNode(config       = config,
                              reference    = prefix["Reference"],
                              input_bam    = input_filename,
                              output_bam   = output_filename,
                              tags         = self.tags,
                              dependencies = prefix["Node"])

        validated_node = IndexAndValidateBAMNode(config, node)
        self.bams["Processed"] = {output_filename : validated_node}


    def _init_reads(self, config, record):
        key = tuple(self.tags[key] for key in ("Target", "SM", "LB", "PU"))
        if key not in _TRIMMED_READS_CACHE:
            _TRIMMED_READS_CACHE[key] = reads.Reads(config, record, record["Options"]["QualityOffset"])
        self.reads = _TRIMMED_READS_CACHE[key]


    def _init_unaligned_lane(self, config, prefix, record):
        aln_key, aln_func = _select_aligner(record["Options"])

        postfix = ["minQ%i" % record["Options"]["Aligners"][aln_key]["MinQuality"]]
        if not record["Options"]["Aligners"][aln_key].get("UseSeed", True):
            postfix.append("noSeed")

        for (key, input_filename) in self.reads.files.iteritems():
            output_filename = os.path.join(self.folder, "%s.%s.bam" % (key.lower(), ".".join(postfix)))
            parameters = {"output_file"  : output_filename,
                          "prefix"       : prefix["Path"],
                          "reference"    : prefix["Reference"],
                          "dependencies" : self.reads.nodes + (prefix["Node"],)}

            alignment_node = aln_func(config, parameters, input_filename, self.tags, record["Options"])
            validated_node = IndexAndValidateBAMNode(config, alignment_node)

            self.bams[key] = {output_filename : validated_node}


def _select_aligner(options):
    key  = options["Aligners"]["Program"]
    if key == "BWA":
        return key, _build_bwa_nodes

    return key, _build_bowtie2_nodes


def _apply_aln_user_parameters(mkfile_params, params, aligners):
    for (param, value) in mkfile_params.iteritems():
        if param.startswith("-"):
            for value in safe_coerce_to_tuple(value):
                for aligner_key in aligners:
                    params.commands[aligner_key].push_parameter(param, value)


def _build_bwa_nodes(config, parameters, input_filename, tags, options):
    if paths.is_paired_end(input_filename):
        input_file_1 = input_file_2 = input_filename
        aln_keys, sam_key = ("aln_1", "aln_2"), "sampe"
    else:
        input_file_1, input_file_2 = input_filename, ""
        aln_keys, sam_key = ("aln",), "samse"

    params = BWANode.customize(input_file_1 = input_file_1.format(Pair = 1),
                               input_file_2 = input_file_2.format(Pair = 2),
                               threads      = config.bwa_max_threads,
                               **parameters)

    for aln_key in aln_keys:
        if not options["Aligners"]["BWA"]["UseSeed"]:
            params.commands[aln_key].set_parameter("-l", 2**16 - 1)
        if options["QualityOffset"] in (64, "Solexa"):
            params.commands[aln_key].set_parameter("-I")
    _apply_aln_user_parameters(options["Aligners"]["BWA"], params, aln_keys)

    pg_tags  = "bwa:CL:%s" % " ".join(map(str, params.commands[aln_keys[0]].call)).replace("%", "%%")
    pg_tags += " | " + " ".join(map(str, params.commands[sam_key].call)).replace("%", "%%")
    params.commands["convert"].push_parameter("--update-pg-tag", pg_tags)

    read_group = "@RG\tID:{ID}\tSM:{SM}\tLB:{LB}\tPU:{PU}\tPL:{PL}\tPG:{PG}".format(**tags)
    params.commands[sam_key].set_parameter("-r", read_group)

    params.commands["filter"].set_parameter('-q', options["Aligners"]["BWA"]["MinQuality"])

    return params.build_node()


def _build_bowtie2_nodes(config, parameters, input_filename, tags, options):
    if paths.is_paired_end(input_filename):
        input_filename_1 = input_filename_2 = input_filename
    else:
        input_filename_1,  input_filename_2 = input_filename, ""

    params = Bowtie2Node.customize(input_file_1    = input_filename_1.format(Pair = 1),
                                   input_file_2    = input_filename_2.format(Pair = 2),
                                   threads         = config.bowtie2_max_threads,
                                   **parameters)

    params.commands["filter"].set_parameter('-q', options["Aligners"]["Bowtie2"]["MinQuality"])
    if options["QualityOffset"] == 64:
        params.commands["aln"].set_parameter("--phred64")
    elif options["QualityOffset"] == 33:
        params.commands["aln"].set_parameter("--phred33")
    else:
        params.commands["aln"].set_parameter("--solexa-quals")
    _apply_aln_user_parameters(options["Aligners"]["Bowtie2"], params, ("aln",))

    pg_tags = "bowtie2:CL:%s" % " ".join(map(str, params.commands["aln"].call)).replace("%", "%%")
    params.commands["convert"].push_parameter("--update-pg-tag", pg_tags)

    params.commands["aln"].set_parameter("--rg-id", tags["ID"])
    for (tag_name, tag_value) in tags.iteritems():
        if tag_name not in ("Target", "ID"):
            params.commands["aln"].push_parameter("--rg", "%s:%s" % (tag_name, tag_value))

    return params.build_node()


def _build_bam_cleanup_nodes(config, target, sample, library, barcode, record):
    tags = {"ID" : library, "SM" : sample, "LB" : library, "PU" : barcode,
            "PL" : record["Options"]["Platform"]}

    results = {}
    for (genome, alignments) in record["Reads"]["BAM"].iteritems():
        results[genome] = {}
        output_dir = os.path.join(config.destination, target, genome, sample, library, barcode)
        for (key, filename) in alignments.iteritems():
            results[genome][key] = {"Node" : node,
                                    "Filename" : output_filename,
                                    "Coverage" : coverages}

    record["Reads"]["BAM"] = results
    return record
