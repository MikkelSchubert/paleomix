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
import itertools
import collections

import pysam

from copy import deepcopy

from pypeline.node import Node, CommandNode, MetaNode
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
from pypeline.atomiccmd.builder import \
     create_customizable_cli_parameters, \
     use_customizable_cli_parameters, \
     AtomicCmdBuilder, \
     apply_options
from pypeline.nodes.samtools import \
     GenotypeNode, \
     TabixIndexNode, \
     FastaIndexNode, \
     MPileupNode
from pypeline.nodes.bedtools import \
    SlopBedNode
from pypeline.nodes.sequences import \
    ExtractReferenceNode


from pypeline.common.fileutils import \
     swap_ext, \
     move_file, \
     add_postfix
from pypeline.common.formats.fasta import FASTA
import pypeline.common.sequences as sequences



class VCFPileupNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, reference, in_bam, in_vcf, outfile, dependencies = ()):
        unicat = AtomicCmdBuilder(["unicat", "%(IN_VCF)s"],
                                  IN_VCF     = in_vcf,
                                  OUT_STDOUT = AtomicCmd.PIPE)

        vcfpileup = AtomicCmdBuilder(["vcf_create_pileup", "%(OUT_PILEUP)s"],
                                     IN_REF       = reference,
                                     IN_BAM       = in_bam,
                                     IN_STDIN     = unicat,
                                     OUT_PILEUP   = outfile,
                                     OUT_TBI      = outfile + ".tbi")
        vcfpileup.add_value("%(IN_BAM)s")
        vcfpileup.set_option("-f", "%(IN_REF)s")

        return {"commands" : {"unicat" : unicat,
                              "pileup" : vcfpileup}}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        commands = [parameters.commands[key].finalize() for key in ("unicat", "pileup")]
        description = "<VCFPileup: '%s' -> '%s'>" % (parameters.in_bam,
                                                     parameters.outfile)
        CommandNode.__init__(self,
                             description  = description,
                             command      = ParallelCmds(commands),
                             dependencies = parameters.dependencies)


class VCFFilterNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, pileup, infile, outfile, regions, dependencies = ()):
        unicat = AtomicCmdBuilder(["unicat", "%(IN_VCF)s"],
                                  IN_VCF     = infile,
                                  OUT_STDOUT = AtomicCmd.PIPE)

        vcffilter = AtomicCmdBuilder(["vcf_filter", "--pileup", "%(IN_PILEUP)s"],
                                     IN_PILEUP = pileup,
                                     IN_STDIN     = unicat,
                                     OUT_STDOUT   = AtomicCmd.PIPE)
        for contig in regions["HomozygousContigs"]:
            vcffilter.add_option("--homozygous-chromosome", contig)

        bgzip = AtomicCmdBuilder(["bgzip"],
                                 IN_STDIN     = vcffilter,
                                 OUT_STDOUT   = outfile)

        return {"commands" : {"unicat" : unicat,
                              "filter" : vcffilter,
                              "bgzip"  : bgzip}}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        commands = [parameters.commands[key].finalize() for key in ("unicat", "filter", "bgzip")]
        description = "<VCFFilter: '%s' -> '%s'>" % (parameters.infile,
                                                     parameters.outfile)
        CommandNode.__init__(self,
                             description  = description,
                             command      = ParallelCmds(commands),
                             dependencies = parameters.dependencies)



class BuildRegionsNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, options, infile, regions, outfile, padding, dependencies = ()):
        params = AtomicCmdBuilder(["bam_genotype_regions", "--padding", padding],
                                  IN_VCFFILE   = infile,
                                  IN_TABIX     = infile + ".tbi",
                                  IN_INTERVALS = regions["BED"],
                                  OUT_STDOUT   = outfile)
        params.set_option("--genotype", "%(IN_VCFFILE)s")
        params.set_option("--intervals", "%(IN_INTERVALS)s")
        if regions["ProteinCoding"]:
            params.set_option("--whole-codon-indels-only")
        if not regions["IncludeIndels"]:
            params.set_option("--ignore-indels")

        return {"command" : params}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()
        description = "<BuildRegions: '%s' -> '%s'>" % (parameters.infile,
                                                        parameters.outfile)
        CommandNode.__init__(self,
                             description  = description,
                             command      = command,
                             dependencies = parameters.dependencies)


class SampleRegionsNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, infile, bedfile, outfile, dependencies = ()):
        params = AtomicCmdBuilder(["bam_sample_regions"],
                                  IN_PILEUP    = infile,
                                  IN_INTERVALS = bedfile,
                                  OUT_STDOUT   = outfile)
        params.set_option("--genotype", "%(IN_PILEUP)s")
        params.set_option("--intervals", "%(IN_INTERVALS)s")

        return {"command" : params}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()

        description = "<SampleRegions: '%s' -> '%s'>" \
            % (parameters.infile, parameters.outfile)

        CommandNode.__init__(self,
                             description  = description,
                             command      = command,
                             dependencies = parameters.dependencies)


_FAI_CACHE = {}
def build_fasta_index_node(reference, dependencies):
    if reference not in _FAI_CACHE:
        _FAI_CACHE[reference] = \
          FastaIndexNode(infile       = reference,
                         dependencies = dependencies)
    return _FAI_CACHE[reference]


_BED_CACHE = {}
def build_regions_nodes(options, regions, padding, dependencies):
    destination = add_postfix(regions["BED"], ".padded_%ibp" % (padding,))

    if not padding:
        return regions["BED"], dependencies

    if destination not in _BED_CACHE:
        faidx_node = build_fasta_index_node(regions["FASTA"], dependencies)
        _BED_CACHE[destination] \
          = SlopBedNode(genome        = regions["FASTA"] + ".fai",
                        infile        = regions["BED"],
                        outfile       = destination,
                        from_start    = padding,
                        from_end      = padding,
                        dependencies  = faidx_node)

    return destination, _BED_CACHE[destination]


def build_genotyping_nodes(options, genotyping, sample, regions, dependencies):
    genotyping = genotyping[regions["Name"]]

    fasta     = regions["Genotypes"][sample["Name"]]
    calls     = swap_ext(fasta, ".vcf.bgz")
    pileups   = swap_ext(fasta, ".vcf.pileup.bgz")
    filtered  = swap_ext(fasta, ".filtered.vcf.bgz")

    padding = genotyping["Padding"]
    infile  = os.path.join(options.samples_root, "%s.%s.bam" % (sample["Name"], regions["Prefix"]))
    if regions["Realigned"]:
        infile = add_postfix(infile, ".realigned")

    slop, node =  build_regions_nodes(options, regions, padding, dependencies)
    genotype = GenotypeNode.customize(reference          = regions["FASTA"],
                                      regions            = slop,
                                      infile             = infile,
                                      outfile            = calls,
                                      dependencies       = node)

    apply_options(genotype.commands["pileup"], genotyping["MPileup"])
    apply_options(genotype.commands["genotype"], genotyping["BCFTools"])
    genotype = genotype.build_node()

    vcfpileup = VCFPileupNode.customize(reference    = regions["FASTA"],
                                        in_bam       = infile,
                                        in_vcf       = calls,
                                        outfile      = pileups,
                                        dependencies = genotype)
    apply_options(vcfpileup.commands["pileup"], genotyping["MPileup"])
    vcfpileup = vcfpileup.build_node()

    vcffilter = VCFFilterNode.customize(infile       = calls,
                                        pileup       = pileups,
                                        outfile      = filtered,
                                        regions      = regions,
                                        dependencies = vcfpileup)

    filter_cfg = genotyping["VCF_Filter"]
    apply_options(vcffilter.commands["filter"], filter_cfg)
    if filter_cfg["MaxReadDepth"]:
        max_depth = filter_cfg["MaxReadDepth"]
        if isinstance(max_depth, dict):
            max_depth = max_depth[sample["Name"]]
        vcffilter.commands["filter"].set_option("--max-read-depth", max_depth)
    vcffilter = vcffilter.build_node()

    tabix    = TabixIndexNode(infile          = filtered,
                              preset          = "vcf",
                              dependencies    = vcffilter)

    builder  = BuildRegionsNode(options      = options,
                                infile       = filtered,
                                regions      = regions,
                                outfile      = fasta,
                                padding      = padding,
                                dependencies = tabix)

    return (builder,)


def build_sampling_nodes(options, genotyping, sample, regions, dependencies):
    genotyping = genotyping[regions["Name"]]

    fasta_file  = regions["Genotypes"][sample["Name"]]
    pileup_file = swap_ext(fasta_file, ".pileup.bgz")

    padding = genotyping["Padding"]
    slop, node =  build_regions_nodes(options, regions, padding, dependencies)

    bam_file = os.path.join(options.samples_root, "%s.%s.bam" % (sample["Name"], regions["Prefix"]))
    if regions["Realigned"]:
        bam_file = add_postfix(bam_file, ".realigned")

    genotype = MPileupNode(reference          = regions["FASTA"],
                           regions            = slop,
                           infile             = bam_file,
                           outfile            = pileup_file,
                           dependencies       = node)
    tabix    = TabixIndexNode(infile          = pileup_file,
                              preset          = "pileup",
                              dependencies    = genotype)

    builder  = SampleRegionsNode(infile       = pileup_file,
                                 bedfile      = regions["BED"],
                                 outfile      = fasta_file,
                                 dependencies = tabix)
    return (builder,)


def build_reference_nodes(options, sample, regions, dependencies):
    output_file = "%s.%s.fasta" % (sample["Name"], regions["Desc"])
    destination = os.path.join(options.destination, "genotypes", output_file)
    faidx_node  = build_fasta_index_node(regions["FASTA"], dependencies)

    node  = ExtractReferenceNode(reference          = regions["FASTA"],
                                 bedfile            = regions["BED"],
                                 outfile            = destination,
                                 dependencies       = faidx_node)
    return (node,)


def build_sample_nodes(options, genotyping, regions_sets, sample, dependencies = ()):
    nodes = []
    for regions in regions_sets.itervalues():
        regions = deepcopy(regions)
        # Enforce homozygous contigs based on gender tag
        regions["HomozygousContigs"] = regions["HomozygousContigs"][sample["Gender"]]

        genotyping_method = sample["GenotypingMethod"].lower()
        if genotyping_method == "reference sequence":
            nodes.extend(build_reference_nodes(options, sample, regions, dependencies))
        elif genotyping_method == "random sampling":
            nodes.extend(build_sampling_nodes(options, genotyping, sample, regions, dependencies))
        elif genotyping_method == "samtools":
            nodes.extend(build_genotyping_nodes(options, genotyping, sample, regions, dependencies))
        else:
            assert False, "Unexpected genotyping method %r for sample %r" % (genotyping_method, sample["Name"])

    return MetaNode(description  = sample["Name"],
                    dependencies = nodes)


def chain(pipeline, options, makefiles):
    destination = options.destination
    for makefile in makefiles:
        regions_sets = makefile["Project"]["Regions"]
        genotyping   = makefile["Genotyping"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        nodes = []
        for sample in makefile["Project"]["Samples"].itervalues():
            nodes.append(build_sample_nodes(options, genotyping, regions_sets, sample, makefile["Nodes"]))

        makefile["Nodes"] = tuple(nodes)
    options.destination = destination
