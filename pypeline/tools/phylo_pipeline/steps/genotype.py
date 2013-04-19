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
import sys
import itertools
import collections

import pysam

from copy import deepcopy

from pypeline.node import Node, CommandNode, MetaNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.nodes.samtools import GenotypeNode, TabixIndexNode, FastaIndexNode, MPileupNode
from pypeline.nodes.bedtools import SlopBedNode


from pypeline.common.fileutils import move_file
import pypeline.common.samwrap as samwrap
import pypeline.common.sequences as sequences
import pypeline.common.formats.fasta as fasta

import common


class BuildRegionsNode(CommandNode):
    def __init__(self, options, settings, infile, interval, outfile, dependencies = ()):
        call = ["bam_genotype_regions",
                "--genotype", "%(IN_VCFFILE)s",
                "--intervals", "%(IN_INTERVALS)s",
                "--padding", settings["Padding"],
                "--min-quality", settings["MinQuality"],
                "--min-read-depth", settings["MinDepth"],
                "--max-read-depth", settings["MaxDepth"],
                "--min-distance-to-indels", settings["MinDistanceToIndels"]]
        if interval["Protein coding"]:
            call.append("--whole-codon-indels-only")
        for contig in interval["Homozygous Contigs"]:
            call.extend(("--homozygous-chromosome", contig))

        prefix = "{Genome}.{Name}".format(**interval)
        intervals = os.path.join(options.intervals_root, prefix + ".bed")
        command = AtomicCmd(call,
                            IN_VCFFILE   = infile,
                            IN_INTERVALS = intervals,
                            OUT_STDOUT   = outfile)

        description = "<BuildRegions: '%s' -> '%s'>" % (infile, outfile)
        CommandNode.__init__(self,
                             description  = description,
                             command      = command,
                             dependencies = dependencies)



class SampleRegionsNode(CommandNode):
    def __init__(self, settings, infile, intervals, outfile, dependencies = ()):
        call = ["bam_sample_regions",
                "--genotype", "%(IN_PILEUP)s",
                "--intervals", "%(IN_INTERVALS)s",
                "--padding", settings["Padding"],
                "--min-distance-to-indels", settings["MinDistanceToIndels"]]

        command = AtomicCmd(call,
                            IN_PILEUP    = infile,
                            IN_INTERVALS = intervals,
                            OUT_STDOUT   = outfile)

        description = "<SampleRegions: '%s' -> '%s'>" \
            % (infile, outfile)

        CommandNode.__init__(self,
                             description  = description,
                             command      = command,
                             dependencies = dependencies)


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
            intervals = samwrap.read_tabix_BED(bedfile).items()
            for (contig, beds) in sorted(intervals):
                beds.sort(key = keyfunc)

                for (gene, gene_beds) in itertools.groupby(beds, lambda x: x.name):
                    gene_beds = tuple(gene_beds)
                    for bed in gene_beds:
                        seqs[(contig, gene)].append(fastafile.fetch(contig, bed.start, bed.end))

                    seq = "".join(seqs[(contig, gene)])
                    if any((bed.strand == "-") for bed in gene_beds):
                        assert all((bed.strand == "-") for bed in gene_beds)
                        seq = sequences.reverse_complement(seq)
                    seqs[(contig, gene)] = seq

        temp_file = os.path.join(temp, "sequences.fasta")
        with open(temp_file, "w") as out_file:
            for ((_, gene), sequence) in sorted(seqs.items()):
                fasta.print_fasta(gene, sequence, out_file)

        move_file(temp_file, self._outfile)





def build_interval_nodes(options, taxa, interval, padding, dependencies):
    prefix = "{Genome}.{Name}".format(**interval)
    source = os.path.join(options.intervals_root, prefix + ".bed")
    genome = os.path.join(options.genomes_root, interval["Genome"] + ".genome")
    destination = os.path.join(options.destination, "genotypes", "%s.%s.bed" % (taxa["Name"], prefix))

    node = SlopBedNode(genome        = genome,
                       infile        = source,
                       outfile       = destination,
                       from_start    = padding,
                       from_end      = padding,
                       dependencies  = dependencies)

    return destination, node


def build_genotyping_nodes(options, genotyping, taxa, interval, dependencies):
    prefix = "{0}.{Genome}.{Name}".format(taxa["Name"], **interval)
    reference = os.path.join(options.genomes_root, interval["Genome"] + ".fasta")
    fasta   = os.path.join(options.destination, "genotypes", prefix + ".fasta")
    pileup  = os.path.join(options.destination, "genotypes", prefix + ".vcf.bgz")

    padding = genotyping["Random"]["Padding"]
    slop, node =  build_interval_nodes(options, taxa, interval, padding, dependencies)
    genotype = GenotypeNode(reference          = reference,
                            regions            = slop,
                            infile             = os.path.join(options.samples_root, "%s.%s.bam" % (taxa["Name"], interval["Genome"])),
                            outfile            = pileup,
                            dependencies       = node)
    tabix    = TabixIndexNode(infile          = pileup,
                              preset          = "vcf",
                              dependencies    = genotype)
    builder  = BuildRegionsNode(options      = options,
                                settings     = genotyping["SAMTools"],
                                infile       = pileup,
                                interval     = interval,
                                outfile      = fasta,
                                dependencies = tabix)
    return (builder,)


def build_sampling_nodes(options, genotyping, taxa, interval, dependencies):
    prefix = "{Genome}.{Name}".format(**interval)
    reference = os.path.join(options.genomes_root, interval["Genome"] + ".fasta")
    destination = os.path.join(options.destination, "genotypes", "%s.%s.fasta" % (taxa["Name"], prefix))
    intervals = os.path.join(options.intervals_root, prefix + ".bed")
    pileup = os.path.join(options.destination, "genotypes", "%s.%s.pileup.bgz" % (taxa["Name"], prefix))

    padding = genotyping["Random"]["Padding"]
    slop, node =  build_interval_nodes(options, taxa, interval, padding, dependencies)
    genotype = MPileupNode(reference          = reference,
                           regions            = slop,
                           infile             = os.path.join(options.samples_root, "%s.%s.bam" % (taxa["Name"], interval["Genome"])),
                           outfile            = pileup,
                           dependencies       = node)
    tabix    = TabixIndexNode(infile          = pileup,
                              preset          = "pileup",
                              dependencies    = genotype)

    builder  = SampleRegionsNode(settings     = genotyping["Random"],
                                 infile       = pileup,
                                 intervals    = intervals,
                                 outfile      = destination,
                                 dependencies = tabix)
    return (builder,)


def build_reference_nodes(options, taxa, interval, dependencies):
    prefix = "{Genome}.{Name}".format(**interval)
    reference = os.path.join(options.genomes_root, interval["Genome"] + ".fasta")
    destination = os.path.join(options.destination, "genotypes", "%s.%s.fasta" % (taxa["Name"], prefix))
    intervals = os.path.join(options.intervals_root, prefix + ".bed")

    faidx = FastaIndexNode(infile               = reference,
                           dependencies         = dependencies)
    node  = ExtractReference(reference          = reference,
                             intervals          = intervals,
                             outfile            = destination,
                             dependencies       = faidx)
    return (node,)


def build_taxa_nodes(options, genotyping, intervals, taxa, dependencies = ()):
    nodes = []
    for interval in intervals.itervalues():
        interval = deepcopy(interval)
        # Override default genome (BAM file) if specified
        interval["Genome"] = common.get_genome_for_interval(interval, taxa)
        # Enforce homozygous contigs based on gender tag
        #        interval[
        interval["Homozygous Contigs"] = interval["Homozygous Contigs"][taxa["Gender"]]

        genotyping_method = taxa.get("Genotyping Method", "samtools").lower()
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
