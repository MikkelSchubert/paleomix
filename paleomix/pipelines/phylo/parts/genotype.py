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

from copy import deepcopy

from paleomix.nodes.samtools import TabixIndexNode, FastaIndexNode, BAMIndexNode
from paleomix.nodes.bedtools import PaddedBedNode
from paleomix.common.fileutils import swap_ext, add_postfix
from paleomix.nodes.commands import (
    VCFFilterNode,
    BuildRegionsNode,
    GenotypeRegionsNode,
)


###############################################################################
###############################################################################

# Caches for nodes shared between multiple tasks
_BAI_CACHE = {}
_FAI_CACHE = {}
_BED_CACHE = {}
_VCF_CACHE = {}


def build_bam_index_node(bamfile):
    """Returns a node generating a BAI index (using SAMTools) for a BAM file;
    the result is cached, to ensure that multiple calls for the same BAM does
    not result in files being clobbered.

    """
    if bamfile not in _BAI_CACHE:
        _BAI_CACHE[bamfile] = BAMIndexNode(infile=bamfile)
    return _BAI_CACHE[bamfile]


def build_fasta_index_node(reference):
    if reference not in _FAI_CACHE:
        _FAI_CACHE[reference] = FastaIndexNode(infile=reference)
    return _FAI_CACHE[reference]


def build_regions_nodes(regions, padding, dependencies=()):
    destination = add_postfix(regions["BED"], ".padded_%ibp" % (padding,))

    if not padding:
        return regions["BED"], dependencies

    if destination not in _BED_CACHE:
        dependencies = list(dependencies)
        dependencies.append(build_fasta_index_node(regions["FASTA"]))
        _BED_CACHE[destination] = PaddedBedNode(
            fai_file=regions["FASTA"] + ".fai",
            infile=regions["BED"],
            outfile=destination,
            amount=padding,
            dependencies=dependencies,
        )

    return destination, (_BED_CACHE[destination],)


def _get_vcf_filter_options(genotyping, sample):
    options = dict(genotyping["VCF_Filter"])
    if options["MaxReadDepth"][sample]:
        options["--max-read-depth"] = options["MaxReadDepth"][sample]
    return options


def build_genotyping_bedfile_nodes(options, genotyping, sample, regions, dependencies):
    bamfile = "%s.%s.bam" % (sample, regions["Prefix"])
    bamfile = os.path.join(options.samples_root, bamfile)
    prefix = regions["Genotypes"][sample]
    padding, bedfile = genotyping["Padding"], None
    if not genotyping["GenotypeEntirePrefix"]:
        bedfile, nodes = build_regions_nodes(regions, padding, dependencies)
        bai_node = build_bam_index_node(bamfile)
        dependencies = nodes + (bai_node,)
    else:
        prefix = os.path.join(
            os.path.dirname(prefix), "%s.%s.TEMP" % (sample, regions["Prefix"])
        )

        dependencies += (build_bam_index_node(bamfile),)

    return prefix, bamfile, bedfile, dependencies


def build_genotyping_nodes_cached(options, genotyping, sample, regions, dependencies):
    """Carries out genotyping, filtering of calls, and indexing of files for a
    given sample and prefix. If the option 'GenotypeEntirePrefix' is enabled,
    the BAM is genotyped once, and each set of RegionsOfInterest simply extract
    the relevant regions during construction of the consensus sequence.

    Parameters:
        options: An options object (c.f. paleomix.pipelines.phylo.config).
        genotyping: Genotyping options defined for a specific set of areas of
                    interest, corresponding to Genotyping:NAME in the makefile.
        sample: The name of the sample to be genotyped.
        egions: A dictionary for a 'RegionsOfInterest' from the makefile.
        dependencies: Depenencies that must be met before genotyping starts.

    Returns a tuple containing the filename of the filtered and tabix-indexed
    VCF file, and the top-level node generating this file. Multiple calls for
    the same BAM and prefix will return the same VCF and nodes if the option
    for 'GenotypeEntirePrefix' is enabled, otherwise each ROI is genotyped
    individiually.

    Output files are generated in ./results/PROJECT/genotyping. If the option
    for 'GenotypeEntirePrefix' is enabled, the following files are generated:
        SAMPLE.PREFIX.vcf.bgz: Unfiltered calls for variant/non-variant sites.
        SAMPLE.PREFIX.filtered.vcf.bgz: Variant calls filtered with vcf_filter.
        SAMPLE.PREFIX.filtered.vcf.bgz.tbi: Tabix index for the filtered VCF.

    If 'GenotypeEntirePrefix' is not enabled for a given ROI, the following
    files are generated for that ROI (see descriptions above):
        SAMPLE.PREFIX.ROI.filtered.vcf.bgz
        SAMPLE.PREFIX.ROI.filtered.vcf.bgz.tbi
        SAMPLE.PREFIX.ROI.vcf.bgz

    In addition, the following files are generated for each set of
    RegionsOfInterest (ROI), regardless of the 'GenotypeEntirePrefix' option:
        SAMPLE.PREFIX.ROI.CDS.fasta: FASTA sequence of each feature in the ROI.
        SAMPLE.PREFIX.ROI.CDS.fasta.fai: FASTA index generated using SAMTools.

    """
    output_prefix, bamfile, bedfile, dependencies = build_genotyping_bedfile_nodes(
        options, genotyping, sample, regions, dependencies
    )

    if (bamfile, output_prefix) in _VCF_CACHE:
        return _VCF_CACHE[(bamfile, output_prefix)]

    calls = swap_ext(output_prefix, ".vcf.bgz")
    filtered = swap_ext(output_prefix, ".filtered.vcf.bgz")

    # 1. Call samtools mpilup | bcftools view on the bam
    genotype = GenotypeRegionsNode(
        reference=regions["FASTA"],
        bedfile=bedfile,
        infile=bamfile,
        outfile=calls,
        mpileup_options=genotyping["MPileup"],
        bcftools_options=genotyping["BCFTools"],
        dependencies=dependencies,
    )

    # 2. Filter all sites using the 'vcf_filter' command
    vcffilter = VCFFilterNode(
        infile=calls,
        outfile=filtered,
        regions=regions,
        options=_get_vcf_filter_options(genotyping, sample),
        dependencies=genotype,
    )

    # 3. Tabix index. This allows random-access to the VCF file when building
    #    the consensus FASTA sequence later in the pipeline.
    tabix = TabixIndexNode(infile=filtered, preset="vcf", dependencies=vcffilter)

    _VCF_CACHE[(bamfile, output_prefix)] = (filtered, tabix)
    return filtered, tabix


def build_genotyping_nodes(options, genotyping, sample, regions, dependencies):
    """Builds the nodes required for genotyping a BAM, in part or in whole.

    By default, only the region of interest (including padding) will be
    genotyped. However, if option 'GenotypeEntirePrefix' is enabled, the entire
    genome is genotyped, and reused between different areas of interest.

    In addition to the files generated by 'build_genotyping_nodes_cached', this
    function generates the following files:
        SAMPLE.PREFIX.ROI.fasta: FASTA containing each named region.
        SAMPLE.PREFIX.ROI.fasta.fai: Index file built using "samtools faidx"

    The function returns a sequence of the top-level nodes generating the files.

    """
    # 1. Get path of the filtered VCF file, and the assosiated node
    filtered, node = build_genotyping_nodes_cached(
        options=options,
        genotyping=genotyping,
        sample=sample,
        regions=regions,
        dependencies=dependencies,
    )

    # 2. Generate consensus sequence from filtered VCF
    output_fasta = regions["Genotypes"][sample]
    builder_options = {}
    if regions["ProteinCoding"]:
        builder_options["--whole-codon-indels-only"] = None
    if not regions["IncludeIndels"]:
        builder_options["--ignore-indels"] = None

    builder = BuildRegionsNode(
        infile=filtered,
        bedfile=regions["BED"],
        outfile=output_fasta,
        padding=genotyping["Padding"],
        options=builder_options,
        dependencies=node,
    )

    # 3. Index sequences to make retrival easier for MSA
    faidx = FastaIndexNode(infile=output_fasta, dependencies=builder)

    return (faidx,)


def build_sample_nodes(options, genotyping, regions_sets, sample, dependencies=()):
    nodes = []
    for regions in regions_sets.values():
        regions = deepcopy(regions)

        # Enforce homozygous contigs based on sex tag
        regions["HomozygousContigs"] = regions["HomozygousContigs"][sample["Sex"]]

        nodes.extend(
            build_genotyping_nodes(
                options=options,
                genotyping=genotyping[regions["Name"]],
                sample=sample["Name"],
                regions=regions,
                dependencies=dependencies,
            )
        )

    return nodes


def chain(pipeline, options, makefiles):
    destination = options.destination
    for makefile in makefiles:
        regions_sets = makefile["Project"]["Regions"]
        genotyping = makefile["Genotyping"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        nodes = []
        for sample in makefile["Project"]["Samples"].values():
            nodes.extend(
                build_sample_nodes(
                    options, genotyping, regions_sets, sample, makefile["Nodes"]
                )
            )

        makefile["Nodes"] = tuple(nodes)
    options.destination = destination
