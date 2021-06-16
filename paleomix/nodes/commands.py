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
import os

import paleomix.tools.bam_stats.coverage as coverage
import paleomix.tools.factory as factory
from paleomix.common.command import AtomicCmd, InputFile, OutputFile, ParallelCmds
from paleomix.common.fileutils import describe_files, move_file, reroot_path
from paleomix.node import CommandNode, Node
from paleomix.nodes.samtools import BCFTOOLS_VERSION, merge_bam_files_command


class CoverageNode(CommandNode):
    def __init__(self, target_name, input_file, output_file, dependencies=()):
        command = factory.new(
            [
                "coverage",
                InputFile(input_file),
                OutputFile(output_file),
                "--target-name",
                target_name,
            ]
        )

        CommandNode.__init__(
            self,
            command=command,
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

    def _run(self, temp):
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
        dependencies=(),
    ):
        command = factory.new(
            [
                "depths",
                InputFile(input_file),
                OutputFile(output_file),
                "--target-name",
                target_name,
            ],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="calculating depth histogram for %s" % (input_file,),
            dependencies=dependencies,
        )


class FilterCollapsedBAMNode(CommandNode):
    def __init__(self, input_bams, output_bam, keep_dupes=True, dependencies=()):
        input_bams = tuple(input_bams)
        if len(input_bams) > 1:
            merge = merge_bam_files_command(input_bams)
            markdup = factory.new(
                "rmdup_collapsed",
                stdin=merge,
                stdout=output_bam,
            )
            command = ParallelCmds([merge, markdup])
        elif len(input_bams) == 1:
            markdup = command = factory.new(
                ["rmdup_collapsed", InputFile(input_bams[0])],
                stdout=output_bam,
            )
        else:
            raise ValueError("empty list of BAM files")

        if not keep_dupes:
            markdup.append("--remove-duplicates")

        CommandNode.__init__(
            self,
            command=command,
            description="filtering merged-read PCR duplicates in %s"
            % (describe_files(input_bams),),
            dependencies=dependencies,
        )


class VCFFilterNode(CommandNode):
    def __init__(self, infile, outfile, regions, options, dependencies=()):
        vcffilter = factory.new(
            ["vcf_filter", InputFile(infile)],
            stdout=AtomicCmd.PIPE,
        )

        vcffilter.merge_options(
            user_options=options,
            fixed_options={
                "--homozygous-chromosome": regions["HomozygousContigs"],
            },
        )

        bgzip = AtomicCmd(["bgzip"], stdin=vcffilter, stdout=outfile)

        CommandNode.__init__(
            self,
            description="filtering VCF records in %s" % (infile,),
            command=ParallelCmds([vcffilter, bgzip]),
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
        mpileup = AtomicCmd(
            ("bcftools", "mpileup", InputFile(infile)),
            stdout=AtomicCmd.PIPE,
            requirements=[BCFTOOLS_VERSION],
        )

        fixed_options = {
            # Ignore read-groups for pileup
            "--ignore-RG": None,
            # Reference sequence (FASTA)
            "--fasta-ref": reference,
            # Output compressed VCF
            "--output-type": "u",
        }

        if bedfile:
            fixed_options["--regions-file"] = InputFile(bedfile)

        mpileup.merge_options(
            user_options=mpileup_options,
            fixed_options=fixed_options,
        )

        genotype = AtomicCmd(
            ("bcftools", "call", "-"),
            stdin=mpileup,
            stdout=outfile,
            requirements=[BCFTOOLS_VERSION],
        )

        genotype.merge_options(
            user_options=bcftools_options,
            fixed_options={
                "--output-type": "z",
            },
        )

        CommandNode.__init__(
            self,
            description="calling genotypes from %s" % (infile,),
            command=ParallelCmds([mpileup, genotype]),
            dependencies=dependencies,
        )


class BuildRegionsNode(CommandNode):
    def __init__(self, infile, bedfile, outfile, padding, options={}, dependencies=()):
        command = factory.new(
            "vcf_to_fasta",
            stdout=outfile,
            extra_files=[
                InputFile(infile + ".tbi"),
            ],
        )

        command.merge_options(
            user_options=options,
            fixed_options={
                "--padding": padding,
                "--genotype": InputFile(infile),
                "--intervals": InputFile(bedfile),
            },
        )

        CommandNode.__init__(
            self,
            description="building FASTA from %s" % (infile,),
            command=command,
            dependencies=dependencies,
        )


class PaddedBedNode(CommandNode):
    def __init__(self, infile, outfile, fai_file, amount=0, dependencies=()):
        command = factory.new(
            [
                ":bedtools",
                "pad",
                "--padding",
                amount,
                InputFile(fai_file),
                InputFile(infile),
            ],
            stdout=outfile,
        )

        CommandNode.__init__(
            self,
            description="padding BED records (%+i) in %s" % (amount, infile),
            command=command,
            dependencies=dependencies,
        )


class FinalizeBAMNode(CommandNode):
    def __init__(
        self,
        in_bams,
        out_passed,
        out_failed,
        out_json,
        options={},
        dependencies=(),
    ):
        in_bams = tuple(in_bams)
        if len(in_bams) > 1:
            merge = merge_bam_files_command(in_bams)
            finalize = factory.new(
                "ngs:finalize_bam",
                stdin=merge,
            )

            command = ParallelCmds([merge, finalize])
        elif len(in_bams) == 1:
            finalize = command = factory.new(
                ["ngs:finalize_bam", InputFile(in_bams[0])]
            )
        else:
            raise ValueError(in_bams)

        finalize.merge_options(
            user_options=options,
            fixed_options={
                "--out-passed": OutputFile(out_passed),
                "--out-failed": OutputFile(out_failed),
                "--out-json": OutputFile(out_json),
            },
        )

        CommandNode.__init__(
            self,
            description="creating finalized BAMs %r and %r"
            % (
                os.path.basename(out_passed),
                os.path.basename(out_failed),
            ),
            command=command,
            dependencies=dependencies,
            threads=options.get("--threads", 1),
        )
