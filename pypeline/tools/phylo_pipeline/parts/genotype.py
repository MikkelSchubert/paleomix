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
from pypeline.nodes.bedtools import SlopBedNode


from pypeline.common.fileutils import \
     move_file, \
     add_postfix
from pypeline.common.formats.fasta import FASTA
import pypeline.common.text as text
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
    def customize(cls, pileup, infile, outfile, interval, dependencies = ()):
        unicat = AtomicCmdBuilder(["unicat", "%(IN_VCF)s"],
                                  IN_VCF     = infile,
                                  OUT_STDOUT = AtomicCmd.PIPE)

        vcffilter = AtomicCmdBuilder(["vcf_filter", "--pileup", "%(IN_PILEUP)s"],
                                     IN_PILEUP = pileup,
                                     IN_STDIN     = unicat,
                                     OUT_STDOUT   = AtomicCmd.PIPE)
        for contig in interval["Homozygous Contigs"]:
            vcffilter.set_option("--homozygous-chromosome", contig)

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
    def customize(cls, options, infile, interval, outfile, padding, dependencies = ()):
        prefix = "{Genome}.{Name}".format(**interval)
        intervals = os.path.join(options.intervals_root, prefix + ".bed")

        params = AtomicCmdBuilder(["bam_genotype_regions"],
                                  IN_VCFFILE   = infile,
                                  IN_TABIX     = infile + ".tbi",
                                  IN_INTERVALS = intervals,
                                  OUT_STDOUT   = outfile)
        params.set_option("--genotype", "%(IN_VCFFILE)s")
        params.set_option("--intervals", "%(IN_INTERVALS)s")
        if interval["Protein coding"]:
            params.set_option("--whole-codon-indels-only")
        if not interval["Include indels"]:
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
    def customize(cls, infile, intervals, outfile, dependencies = ()):
        params = AtomicCmdBuilder(["bam_sample_regions"],
                                  IN_PILEUP    = infile,
                                  IN_INTERVALS = intervals,
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


class ExtractReference(Node):
    def __init__(self, reference, intervals, outfile, dependencies = ()):
        self._reference = reference
        self._intervals = intervals
        self._outfile   = outfile

        description = "<ExtractReference: '%s' -> '%s'>" \
            % (reference, outfile)
        Node.__init__(self,
                      description  = description,
                      input_files  = [intervals],
                      output_files = [outfile],
                      dependencies = dependencies)


    def _run(self, _config, temp):
        def keyfunc(bed):
            return (bed.contig, bed.name, bed.start)

        fastafile = pysam.Fastafile(self._reference)
        seqs = collections.defaultdict(list)
        with open(self._intervals) as bedfile:
            intervals = text.parse_lines_by_contig(bedfile, pysam.asBed()).items()
            for (contig, beds) in sorted(intervals):
                beds.sort(key = keyfunc)

                for (gene, gene_beds) in itertools.groupby(beds, lambda x: x.name):
                    gene_beds = tuple(gene_beds)
                    for bed in gene_beds:
                        # TODO: Check that region was fetched!
                        seqs[(contig, gene)].append(fastafile.fetch(contig, bed.start, bed.end))

                    seq = "".join(seqs[(contig, gene)])
                    if any((bed.strand == "-") for bed in gene_beds):
                        assert all((bed.strand == "-") for bed in gene_beds)
                        seq = sequences.reverse_complement(seq)
                    seqs[(contig, gene)] = seq

        temp_file = os.path.join(temp, "sequences.fasta")
        with open(temp_file, "w") as out_file:
            for ((_, gene), sequence) in sorted(seqs.items()):
                FASTA(gene, None, sequence).write(out_file)

        move_file(temp_file, self._outfile)


_FAI_CACHE = {}
def build_fasta_index_node(reference, dependencies):
    if reference not in _FAI_CACHE:
        _FAI_CACHE[reference] = \
          FastaIndexNode(infile       = reference,
                         dependencies = dependencies)
    return _FAI_CACHE[reference]


_BED_CACHE = {}
def build_interval_nodes(options, interval, padding, dependencies):
    prefix = "{Genome}.{Name}".format(**interval)
    source = os.path.join(options.intervals_root, prefix + ".bed")
    reference   = os.path.join(options.genomes_root, interval["Genome"] + ".fasta")
    destination = add_postfix(source, ".padded_%ibp" % (padding,))

    if not padding:
        return source, dependencies

    if destination not in _BED_CACHE:
        faidx_node = build_fasta_index_node(reference, dependencies)
        _BED_CACHE[destination] \
          = SlopBedNode(genome        = reference + ".fai",
                        infile        = source,
                        outfile       = destination,
                        from_start    = padding,
                        from_end      = padding,
                        dependencies  = faidx_node)

    return destination, _BED_CACHE[destination]


def build_genotyping_nodes(options, genotyping, taxa, interval, dependencies):
    prefix = "{0}.{Genome}.{Name}".format(taxa["Name"], **interval)
    reference = os.path.join(options.genomes_root, interval["Genome"] + ".fasta")
    fasta     = os.path.join(options.destination, "genotypes", prefix + ".fasta")
    calls     = os.path.join(options.destination, "genotypes", prefix + ".vcf.bgz")
    pileups   = os.path.join(options.destination, "genotypes", prefix + ".vcf.pileup.bgz")
    filtered  = os.path.join(options.destination, "genotypes", prefix + ".filtered.vcf.bgz")

    padding = genotyping["Padding"]
    infile  = os.path.join(options.samples_root, "%s.%s.bam" % (taxa["Name"], interval["Genome"]))
    slop, node =  build_interval_nodes(options, interval, padding, dependencies)
    genotype = GenotypeNode.customize(reference          = reference,
                                      regions            = slop,
                                      infile             = infile,
                                      outfile            = calls,
                                      dependencies       = node)

    apply_options(genotype.commands["pileup"], genotyping["MPileup"])
    apply_options(genotype.commands["genotype"], genotyping["BCFTools"])
    genotype = genotype.build_node()

    vcfpileup = VCFPileupNode.customize(reference    = reference,
                                        in_bam       = infile,
                                        in_vcf       = calls,
                                        outfile      = pileups,
                                        dependencies = genotype)
    apply_options(vcfpileup.commands["pileup"], genotyping["MPileup"])
    vcfpileup = vcfpileup.build_node()

    vcffilter = VCFFilterNode.customize(infile       = calls,
                                        pileup       = pileups,
                                        outfile      = filtered,
                                        interval     = interval,
                                        dependencies = vcfpileup)

    filter_cfg = genotyping["VCF_Filter"]
    apply_options(vcffilter.commands["filter"], filter_cfg)
    if filter_cfg["MaxReadDepth"]:
        max_depth = filter_cfg["MaxReadDepth"]
        if isinstance(max_depth, dict):
            max_depth = max_depth[taxa["Name"]]
        vcffilter.commands["filter"].set_option("--max-read-depth", max_depth)
    vcffilter = vcffilter.build_node()

    tabix    = TabixIndexNode(infile          = filtered,
                              preset          = "vcf",
                              dependencies    = vcffilter)

    builder  = BuildRegionsNode(options      = options,
                                infile       = filtered,
                                interval     = interval,
                                outfile      = fasta,
                                padding      = padding,
                                dependencies = tabix)

    return (builder,)


def build_sampling_nodes(options, genotyping, taxa, interval, dependencies):
    prefix = "{Genome}.{Name}".format(**interval)
    reference = os.path.join(options.genomes_root, interval["Genome"] + ".fasta")
    destination = os.path.join(options.destination, "genotypes", "%s.%s.fasta" % (taxa["Name"], prefix))
    intervals = os.path.join(options.intervals_root, prefix + ".bed")
    pileup = os.path.join(options.destination, "genotypes", "%s.%s.pileup.bgz" % (taxa["Name"], prefix))

    padding = genotyping["Padding"]
    slop, node =  build_interval_nodes(options, interval, padding, dependencies)
    genotype = MPileupNode(reference          = reference,
                           regions            = slop,
                           infile             = os.path.join(options.samples_root, "%s.%s.bam" % (taxa["Name"], interval["Genome"])),
                           outfile            = pileup,
                           dependencies       = node)
    tabix    = TabixIndexNode(infile          = pileup,
                              preset          = "pileup",
                              dependencies    = genotype)

    builder  = SampleRegionsNode(infile       = pileup,
                                 intervals    = intervals,
                                 outfile      = destination,
                                 dependencies = tabix)
    return (builder,)


def build_reference_nodes(options, taxa, interval, dependencies):
    prefix = "{Genome}.{Name}".format(**interval)
    reference = os.path.join(options.genomes_root, taxa["Name"] + ".fasta")
    destination = os.path.join(options.destination, "genotypes", "%s.%s.fasta" % (taxa["Name"], prefix))
    intervals = os.path.join(options.intervals_root, prefix + ".bed")
    faidx_node = build_fasta_index_node(reference, dependencies)

    node  = ExtractReference(reference          = reference,
                             intervals          = intervals,
                             outfile            = destination,
                             dependencies       = faidx_node)
    return (node,)


def build_taxa_nodes(options, genotyping, intervals, taxa, dependencies = ()):
    nodes = []
    for interval in intervals.itervalues():
        interval = deepcopy(interval)
        # Enforce homozygous contigs based on gender tag
        interval["Homozygous Contigs"] = interval["Homozygous Contigs"][taxa["Gender"]]

        genotyping_method = taxa["Genotyping Method"]
        if genotyping_method == "reference sequence":
            nodes.extend(build_reference_nodes(options, taxa, interval, dependencies))
        elif genotyping_method == "random sampling":
            nodes.extend(build_sampling_nodes(options, genotyping, taxa, interval, dependencies))
        elif genotyping_method == "samtools":
            nodes.extend(build_genotyping_nodes(options, genotyping, taxa, interval, dependencies))

    return MetaNode(description  = taxa["Name"],
                    dependencies = nodes)


def chain(pipeline, options, makefiles):
    destination = options.destination
    for makefile in makefiles:
        intervals  = makefile["Project"]["Intervals"]
        genotyping = makefile["Genotyping"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        nodes = []
        for taxa in makefile["Project"]["Taxa"].itervalues():
            nodes.append(build_taxa_nodes(options, genotyping, intervals, taxa, makefile["Nodes"]))

        makefile["Nodes"] = tuple(nodes)
    options.destination = destination
