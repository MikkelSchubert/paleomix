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
"""Implements nodes for calling PALEOMIX commands.

Each node is equivalent to a particular command:
    $ paleomix [...]
"""
from paleomix.node import CommandNode, Node
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import ParallelCmds
from paleomix.atomiccmd.builder import (
    AtomicCmdBuilder,
    apply_options,
)
from paleomix.common.fileutils import describe_files, reroot_path, move_file
from paleomix.nodes.samtools import merge_bam_files_command, BCFTOOLS_VERSION

import paleomix.tools.bam_stats.coverage as coverage
import paleomix.tools.factory as factory


class CoverageNode(CommandNode):
    def __init__(
        self, target_name, input_file, output_file, regions_file=None, dependencies=()
    ):
        builder = factory.new("coverage")
        builder.add_value("%(IN_BAM)s")
        builder.add_value("%(OUT_FILE)s")
        builder.set_option("--target-name", target_name)
        builder.set_kwargs(IN_BAM=input_file, OUT_FILE=output_file)

        if regions_file:
            builder.set_option("--regions-file", "%(IN_REGIONS)s")
            builder.set_kwargs(IN_REGIONS=regions_file)

        CommandNode.__init__(
            self,
            command=builder.finalize(),
            description="calculating coverage for %s" % (input_file,),
            dependencies=dependencies,
        )


class MergeCoverageNode(Node):
    def __init__(self, input_files, output_file, dependencies=()):
        self._output_file = output_file

        Node.__init__(
            self,
            description="merging coverage tables from %s"
            % (describe_files(input_files),),
            input_files=input_files,
            output_files=self._output_file,
            dependencies=dependencies,
        )

    def _run(self, _config, temp):
        table = {}
        for filename in self.input_files:
            coverage.read_table(table, filename)

        coverage.write_table(table, reroot_path(temp, self._output_file))
        move_file(reroot_path(temp, self._output_file), self._output_file)


class DepthHistogramNode(CommandNode):
    def __init__(
        self,
        target_name,
        input_file,
        output_file,
        prefix,
        regions_file=None,
        dependencies=(),
    ):
        index_format = regions_file and prefix["IndexFormat"]

        builder = factory.new("depths")
        builder.add_value("%(IN_BAM)s")
        builder.add_value("%(OUT_FILE)s")
        builder.set_option("--target-name", target_name)
        builder.set_kwargs(OUT_FILE=output_file, IN_BAM=input_file)

        if regions_file:
            builder.set_option("--regions-file", "%(IN_REGIONS)s")
            builder.set_kwargs(
                IN_REGIONS=regions_file, TEMP_IN_INDEX=input_file + index_format
            )

        CommandNode.__init__(
            self,
            command=builder.finalize(),
            description="calculating depth histogram for %s" % (input_file,),
            dependencies=dependencies,
        )


class FilterCollapsedBAMNode(CommandNode):
    def __init__(
        self, config, input_bams, output_bam, keep_dupes=True, dependencies=()
    ):
        merge = merge_bam_files_command(input_bams)

        builder = factory.new("rmdup_collapsed")
        builder.set_kwargs(IN_STDIN=merge, OUT_STDOUT=output_bam)

        if not keep_dupes:
            builder.set_option("--remove-duplicates")

        CommandNode.__init__(
            self,
            command=ParallelCmds([merge, builder.finalize()]),
            description="filtering collapsed PCR duplicates in %s"
            % (describe_files(merge.input_files),),
            dependencies=dependencies,
        )


class VCFFilterNode(CommandNode):
    def __init__(self, infile, outfile, regions, options, dependencies=()):
        vcffilter = factory.new("vcf_filter")
        vcffilter.add_value("%(IN_VCF)s")

        for contig in regions["HomozygousContigs"]:
            vcffilter.add_option("--homozygous-chromosome", contig)
        vcffilter.set_kwargs(IN_VCF=infile, OUT_STDOUT=AtomicCmd.PIPE)

        apply_options(vcffilter, options)

        bgzip = AtomicCmdBuilder(["bgzip"], IN_STDIN=vcffilter, OUT_STDOUT=outfile)

        CommandNode.__init__(
            self,
            description="filtering VCF records in %s" % (infile,),
            command=ParallelCmds([vcffilter.finalize(), bgzip.finalize()]),
            dependencies=dependencies,
        )


class GenotypeRegionsNode(CommandNode):
    def __init__(
        self,
        reference,
        infile,
        bedfile,
        outfile,
        mpileup_options={},
        bcftools_options={},
        dependencies=(),
    ):
        mpileup = AtomicCmdBuilder(
            ("bcftools", "mpileup", "%(IN_BAMFILE)s"),
            IN_BAMFILE=infile,
            IN_INTERVALS=bedfile,
            OUT_STDOUT=AtomicCmd.PIPE,
            CHECK_VERSION=BCFTOOLS_VERSION,
        )

        # Ignore read-groups for pileup
        mpileup.add_option("--ignore-RG")
        # Reference sequence (FASTA)
        mpileup.add_option("--fasta-ref", reference)
        # Output compressed VCF
        mpileup.add_option("--output-type", "u")

        if bedfile:
            mpileup.set_option("--regions-file", "%(IN_INTERVALS)s")

        apply_options(mpileup, mpileup_options)

        genotype = AtomicCmdBuilder(
            ("bcftools", "call", "-"),
            IN_STDIN=mpileup,
            IN_BAMFILE=infile,
            OUT_STDOUT=outfile,
            CHECK_VERSION=BCFTOOLS_VERSION,
        )

        genotype.set_option("--output-type", "z")

        apply_options(genotype, bcftools_options)

        CommandNode.__init__(
            self,
            description="calling genotypes from %s" % (infile,),
            command=ParallelCmds([mpileup.finalize(), genotype.finalize()]),
            dependencies=dependencies,
        )


class BuildRegionsNode(CommandNode):
    def __init__(self, infile, bedfile, outfile, padding, options={}, dependencies=()):
        params = factory.new("vcf_to_fasta")
        params.set_option("--padding", padding)
        params.set_option("--genotype", "%(IN_VCFFILE)s")
        params.set_option("--intervals", "%(IN_INTERVALS)s")

        params.set_kwargs(
            IN_VCFFILE=infile,
            IN_TABIX=infile + ".tbi",
            IN_INTERVALS=bedfile,
            OUT_STDOUT=outfile,
        )

        apply_options(params, options)

        CommandNode.__init__(
            self,
            description="building FASTA from %s" % (infile,),
            command=params.finalize(),
            dependencies=dependencies,
        )


def _apply_samtools_options(builder, options, argument):
    for (key, value) in dict(options).items():
        sam_argument = key
        if value is not None:
            sam_argument = "%s=%s" % (key, value)

        builder.add_option(argument, sam_argument, sep="=")
