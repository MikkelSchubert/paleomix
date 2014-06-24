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


import pypeline.tools.bam_pipeline.paths as paths

from pypeline.common.utilities import \
    safe_coerce_to_tuple
from pypeline.atomiccmd.builder import \
    apply_options
from pypeline.nodes.bwa import \
    BWANode
from pypeline.nodes.bowtie2 import \
    Bowtie2Node
from pypeline.tools.bam_pipeline.parts import \
    Reads
from pypeline.tools.bam_pipeline.nodes import \
    CleanupBAMNode, \
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
        self.tags["PU"] = self.tags["PU_src"]
        self.tags["PG"] = self.tags["PG"].lower()
        self.folder  = os.path.join(config.destination, tags["Target"], prefix["Name"],
                                    tags["SM"], tags["LB"], tags["PU_cur"])

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

        validated_node = IndexAndValidateBAMNode(config, prefix, node)
        self.bams["Processed"] = {output_filename : validated_node}


    def _init_reads(self, config, record):
        key = tuple(self.tags[key] for key in ("Target", "SM", "LB", "PU_cur"))
        if key not in _TRIMMED_READS_CACHE:
            _TRIMMED_READS_CACHE[key] = Reads(config, record, record["Options"]["QualityOffset"])
        self.reads = _TRIMMED_READS_CACHE[key]


    def _init_unaligned_lane(self, config, prefix, record):
        aln_key, aln_func = _select_aligner(record["Options"])
        prefix_key = "Node:%s" % (aln_key,)

        postfix = ["minQ%i" % record["Options"]["Aligners"][aln_key]["MinQuality"]]
        if not record["Options"]["Aligners"][aln_key].get("UseSeed", True):
            postfix.append("noSeed")

        for (key, input_filename) in self.reads.files.iteritems():
            # Common parameters between BWA / Bowtie2
            output_filename = os.path.join(self.folder, "%s.%s.bam" % (key.lower(), ".".join(postfix)))
            parameters = {"output_file"  : output_filename,
                          "prefix"       : prefix["Path"],
                          "reference"    : prefix["Reference"],
                          "input_file_1" : input_filename,
                          "input_file_2" : None,
                          "dependencies" : self.reads.nodes + (prefix[prefix_key],)}

            if paths.is_paired_end(input_filename):
                parameters["input_file_1"] = input_filename.format(Pair = 1)
                parameters["input_file_2"] = input_filename.format(Pair = 2)

            alignment_obj  = aln_func(config           = config,
                                      parameters       = parameters,
                                      tags             = self.tags,
                                      options          = record["Options"])

            alignment_opts = record["Options"]["Aligners"][aln_key]
            alignment_obj.commands["convert"].set_option('-q', alignment_opts["MinQuality"])
            if alignment_opts["FilterUnmappedReads"]:
                alignment_obj.commands["convert"].set_option('-F', "0x4")

            alignment_node = alignment_obj.build_node()
            validated_node = IndexAndValidateBAMNode(config, prefix, alignment_node)

            self.bams[key] = {output_filename : validated_node}


def _select_aligner(options):
    key  = options["Aligners"]["Program"]
    if key == "BWA":
        return key, _bwa_build_nodes
    elif key == "Bowtie2":
        return key, _bowtie2_build_nodes
    assert False


def _build_mapper_cl_tag(options, cli_tag):
    samtools = ["|", "samtools", "view", "-q", options["MinQuality"]]
    if options["FilterUnmappedReads"]:
        samtools.extend(("-F", "0x4"))
    samtools.append("...")

    cli_tag.extend(samtools)
    cli_tag.extend(("|", "samtools", "fixmate", "..."))
    cli_tag.extend(("|", "samtools", "calmd",   "..."))


def _set_rg_tags(command, rg_tags, pg_tags):
    command.add_option("--update-pg-tag", pg_tags)
    command.set_option("--rg-id", rg_tags["ID"])
    for tag_name in ("SM", "LB", "PU", "PL", "PG"):
        tag_value = "%s:%s" % (tag_name, rg_tags[tag_name])
        command.add_option("--rg", tag_value)


class ParamCollector:
    def __init__(self, prefix = (), postfix = ()):
        self._prefix  = safe_coerce_to_tuple(prefix)
        self._postfix = safe_coerce_to_tuple(postfix)
        self._added = []
        self._set   = {}

    def add_option(self, key, value):
        self._added.extend((key, value))

    def set_option(self, key, value = None):
        self._set[key] = value

    def pop_option(self, _key):
        pass

    def get_result(self):
        result = list(self._prefix)
        result.extend(self._added)
        for (key, value) in sorted(self._set.items()):
            result.append(key)
            if value is not None:
                result.append(value)
        result.extend(self._postfix)
        return result


################################################################################
################################################################################
## BWA

def _bwa_aln_parameters(options):
    if not options["Aligners"]["BWA"]["UseSeed"]:
        yield ("-l", 2**16 - 1)

    if options["QualityOffset"] in (64, "Solexa"):
        yield ("-I", None)

    for (key, value) in options["Aligners"]["BWA"].iteritems():
        yield (key, value)


def _bwa_build_cl_tag(options):
    # Build summary of parameters used by alignment, only including
    # parameters that affect the output of BWA (as far as possible)
    algorithm = options["Aligners"]["BWA"]["Algorithm"].lower()
    algorithm = "aln" if algorithm == "backtrack" else algorithm

    cli_tag = ParamCollector(("bwa", algorithm), "...")
    apply_options(cli_tag, _bwa_aln_parameters(options))
    cli_tag = cli_tag.get_result()

    _build_mapper_cl_tag(options["Aligners"]["BWA"], cli_tag)

    return " ".join(map(str, cli_tag)).replace("%", "%%")


def _bwa_build_nodes(config, parameters, tags, options):
    algorithm = options["Aligners"]["BWA"]["Algorithm"].lower()

    params = BWANode(threads=config.bwa_max_threads,
                     algorithm=algorithm,
                     **parameters)

    parameters = dict(_bwa_aln_parameters(options))
    # "aln" is used by SE backtrack, mem, and sw; _1 and _2 by PE backtrack
    for aln_key in ("aln", "aln_1", "aln_2"):
        if aln_key in params.commands:
            apply_options(params.commands[aln_key], parameters)

    pg_tags = "bwa:CL:%s" % (_bwa_build_cl_tag(options),)
    _set_rg_tags(params.commands["convert"], tags, pg_tags)

    return params


###############################################################################
###############################################################################
## Bowtie2:

def _bowtie2_aln_parameters(options):
    if options["QualityOffset"] == 64:
        yield ("--phred64", None)
    elif options["QualityOffset"] == 33:
        yield ("--phred33", None)
    else:
        yield ("--solexa-quals", None)

    for (key, value) in options["Aligners"]["Bowtie2"].iteritems():
        if key.startswith("-"):
            yield (key, value)


def _bowtie2_build_cl_tag(options):
    # Build summary of parameters used by alignment, only including
    # parameters that affect the output of bowtie2 (as far as possible)
    cli_tag = ParamCollector("bowtie2", "...")
    apply_options(cli_tag, _bowtie2_aln_parameters(options))
    cli_tag = cli_tag.get_result()

    _build_mapper_cl_tag(options["Aligners"]["Bowtie2"], cli_tag)

    return " ".join(map(str, cli_tag)).replace("%", "%%")


def _bowtie2_build_nodes(config, parameters, tags, options):
    params = Bowtie2Node.customize(threads         = config.bowtie2_max_threads,
                                   **parameters)

    apply_options(params.commands["aln"], _bowtie2_aln_parameters(options))

    pg_tags = "bowtie2:CL:%s" % (_bowtie2_build_cl_tag(options),)
    _set_rg_tags(params.commands["convert"], tags, pg_tags)

    return params
