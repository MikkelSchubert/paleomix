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
"""Implements nodes for calling PALEOMIX commands.

Each node is equivalent to a particular command:
    $ paleomix [...]
"""
import os
import gzip

from paleomix.node import \
    CommandNode, \
    Node
from paleomix.atomiccmd.command import \
    AtomicCmd
from paleomix.atomiccmd.sets import \
    ParallelCmds
from paleomix.common.fileutils import \
    describe_files
from paleomix.nodes.picard import \
    MultiBAMInput, \
    MultiBAMInputNode
from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder, \
    create_customizable_cli_parameters, \
    use_customizable_cli_parameters
from paleomix.common.fileutils import \
    reroot_path, \
    move_file

from paleomix.nodes.samtools import \
    SAMTOOLS_VERSION, \
    SAMTOOLS_VERSION_0119, \
    BCFTOOLS_VERSION_0119

import paleomix.tools.bam_stats.coverage as coverage
import paleomix.tools.factory as factory


class DuplicateHistogramNode(MultiBAMInputNode):
    """Node for calling the 'paleomix duphist' command.

    Takes 1 or more BAMs as imput, requiring a config object in order to run
    MergeSamFiles.jar to merge these files. The output is a histogram of PCR
    duplicate counts, usable as input for the 'preseq' tool.
    """

    def __init__(self, config, input_files, output_file, dependencies=()):
        bam_input = MultiBAMInput(config, input_files, indexed=False)
        duphist_command = factory.new("duphist")
        duphist_command.add_value('%(TEMP_IN_BAM)s')
        duphist_command.set_kwargs(OUT_STDOUT=output_file)
        bam_input.setup(duphist_command)
        duphist_command = duphist_command.finalize()

        commands = ParallelCmds(bam_input.commands + [duphist_command])

        description = "<DuplicateHistogram: %s -> %r>" \
            % (describe_files(input_files), output_file)
        MultiBAMInputNode.__init__(self,
                                   bam_input=bam_input,
                                   command=commands,
                                   description=description,
                                   dependencies=dependencies)


class CoverageNode(CommandNode):
    def __init__(self, config, target_name, input_file, output_file,
                 regions_file=None, dependencies=()):
        builder = factory.new("coverage")
        builder.add_value("%(IN_BAM)s")
        builder.add_value("%(OUT_FILE)s")
        builder.set_option("--target-name", target_name)
        builder.set_kwargs(IN_BAM=input_file,
                           OUT_FILE=output_file)

        if regions_file:
            builder.set_option('--regions-file', '%(IN_REGIONS)s')
            builder.set_kwargs(IN_REGIONS=regions_file)

        description = "<Coverage: %s -> '%s'>" % (input_file, output_file)
        CommandNode.__init__(self,
                             command=builder.finalize(),
                             description=description,
                             dependencies=dependencies)


class MergeCoverageNode(Node):
    def __init__(self, input_files, output_file, dependencies=()):
        self._output_file = output_file

        Node.__init__(self,
                      description="<MergeCoverage: %s -> '%s'>"
                      % (describe_files(input_files), self._output_file),
                      input_files=input_files,
                      output_files=self._output_file,
                      dependencies=dependencies)

    def _run(self, _config, temp):
        table = {}
        for filename in self.input_files:
            coverage.read_table(table, filename)

        coverage.write_table(table, reroot_path(temp, self._output_file))
        move_file(reroot_path(temp, self._output_file), self._output_file)


class DepthHistogramNode(MultiBAMInputNode):
    def __init__(self, config, target_name, input_files, output_file,
                 regions_file=None, dependencies=()):
        bam_input = MultiBAMInput(config, input_files,
                                  indexed=bool(regions_file))
        if len(bam_input.files) > 1 and regions_file:
            raise ValueError("DepthHistogram for regions require single, "
                             "indexed input BAM file.")

        builder = factory.new("depths")
        builder.add_value("%(TEMP_IN_BAM)s")
        builder.add_value("%(OUT_FILE)s")
        builder.set_option("--target-name", target_name)
        builder.set_kwargs(OUT_FILE=output_file)
        bam_input.setup(builder)

        if regions_file:
            builder.set_option('--regions-file', '%(IN_REGIONS)s')
            builder.set_kwargs(IN_REGIONS=regions_file)

        command = ParallelCmds(bam_input.commands + [builder.finalize()])
        description = "<DepthHistogram: %s -> '%s'>" \
            % (describe_files(bam_input.files), output_file)

        MultiBAMInputNode.__init__(self,
                                   bam_input=bam_input,
                                   command=command,
                                   description=description,
                                   dependencies=dependencies)


class FilterCollapsedBAMNode(MultiBAMInputNode):
    def __init__(self, config, input_bams, output_bam, keep_dupes=True,
                 dependencies=()):
        bam_input = MultiBAMInput(config, input_bams, indexed=False)

        builder = factory.new("rmdup_collapsed")
        builder.add_value("%(TEMP_IN_BAM)s")
        builder.set_kwargs(OUT_STDOUT=output_bam)
        bam_input.setup(builder)

        if not keep_dupes:
            builder.set_option("--remove-duplicates")

        filteruniq = builder.finalize()
        command = ParallelCmds(bam_input.commands + [filteruniq])
        description = "<FilterCollapsedBAM: %s>" \
            % (describe_files(bam_input.files),)
        MultiBAMInputNode.__init__(self,
                                   bam_input=bam_input,
                                   command=command,
                                   description=description,
                                   dependencies=dependencies)


class VCFPileupNode(CommandNode):
    """Collects heterozygous SNPs from a VCF file, and generates a bgzipped
    pileup for those sites containing the SNPs.

    The resulting pileup is read by 'paleomix vcf_filter'; this allows
    filtering based on the frequency of the minority SNP, since this is not
    reported in the VCF.
    """

    @create_customizable_cli_parameters
    def customize(cls, reference, infile_bam, infile_vcf, outfile,
                  dependencies=()):
        params = factory.new("genotype")
        params.add_value("%(IN_BAMFILE)s")
        params.add_value("%(OUT_PILEUP)s")
        params.set_option("--bedfile", "%(TEMP_IN_INTERVALS)s")
        params.set_option("--pileup-only")
        # Ignore read-groups for pileup
        params.add_option("--mpileup-argument", "-R", sep="=")
        # Reference sequence (FASTA)
        params.add_option("--mpileup-argument",
                          "-f=%s" % (reference,), sep="=")

        params.set_kwargs(IN_BAMFILE=infile_bam,
                          TEMP_IN_INTERVALS="heterozygous_snps.bed",
                          # Automatically remove this file
                          TEMP_OUT_INTERVALS="heterozygous_snps.bed",
                          OUT_PILEUP=outfile,
                          CHECK_SAMTOOLS=SAMTOOLS_VERSION)

        return {"command": params}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._in_vcf = parameters.infile_vcf
        command = parameters.command.finalize()
        description = "<VCFPileup: '%s' -> '%s'>" \
            % (parameters.infile_vcf,
               parameters.outfile)

        CommandNode.__init__(self,
                             description=description,
                             command=command,
                             dependencies=parameters.dependencies)

    def _run(self, config, temp):
        with gzip.open(self._in_vcf) as handle:
            with open(os.path.join(temp, "heterozygous_snps.bed"), "w") as bed:
                for line in handle:
                    if line.startswith("#"):
                        continue

                    fields = line.split("\t", 5)
                    if "," in fields[4]:
                        pos = int(fields[1])
                        bed.write("%s\t%i\t%i\n" % (fields[0], pos - 1, pos))

        CommandNode._run(self, config, temp)


class VCFFilterNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, pileup, infile, outfile, regions, dependencies=()):
        cat = factory.new("cat")
        cat.add_value("%(IN_VCF)s")
        cat.set_kwargs(IN_VCF=infile,
                       OUT_STDOUT=AtomicCmd.PIPE)

        vcffilter = factory.new("vcf_filter")
        vcffilter.add_option("--pileup", "%(IN_PILEUP)s")
        for contig in regions["HomozygousContigs"]:
            vcffilter.add_option("--homozygous-chromosome", contig)
        vcffilter.set_kwargs(IN_PILEUP=pileup,
                             IN_STDIN=cat,
                             OUT_STDOUT=AtomicCmd.PIPE)

        bgzip = AtomicCmdBuilder(["bgzip"],
                                 IN_STDIN=vcffilter,
                                 OUT_STDOUT=outfile)

        return {"commands": {"cat": cat,
                             "filter": vcffilter,
                             "bgzip": bgzip}}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        commands = [parameters.commands[key].finalize()
                    for key in ("cat", "filter", "bgzip")]

        description = "<VCFFilter: '%s' -> '%s'>" % (parameters.infile,
                                                     parameters.outfile)
        CommandNode.__init__(self,
                             description=description,
                             command=ParallelCmds(commands),
                             dependencies=parameters.dependencies)


class GenotypeRegionsNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, reference, infile, bedfile, outfile,
                  pileup_only=False, nbatches=1, dependencies=()):
        params = factory.new("genotype")
        params.add_value("%(IN_BAMFILE)s")
        params.add_value("%(OUT_VCFFILE)s")
        params.set_option("--nbatches", nbatches)

        if bedfile:
            params.set_option("--bedfile", "%(IN_INTERVALS)s")

        if pileup_only:
            params.set_option("--pileup-only")
            # Ignore read-groups for pileup
            params.add_option("--mpileup-argument", "-R", sep="=")

        # Reference sequence (FASTA)
        params.add_option("--mpileup-argument",
                          "-f=%s" % (reference,), sep="=")

        params.set_kwargs(IN_BAMFILE=infile,
                          IN_INTERVALS=bedfile,
                          OUT_VCFFILE=outfile,
                          CHECK_SAMTOOLS=SAMTOOLS_VERSION_0119,
                          CHECK_BCFTOOLS=BCFTOOLS_VERSION_0119)

        return {"command": params}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()
        invokation = " (%s%i thread(s))" \
            % ("pileup; " if parameters.pileup_only else "",
               parameters.nbatches)
        description = "<GenotypeRegions%s: '%s' -> '%s'>" \
            % (invokation,
               parameters.infile,
               parameters.outfile)

        CommandNode.__init__(self,
                             description=description,
                             command=command,
                             threads=parameters.nbatches,
                             dependencies=parameters.dependencies)


class BuildRegionsNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, infile, bedfile, outfile, padding, dependencies=()):
        params = factory.new("vcf_to_fasta")
        params.set_option("--padding", padding)
        params.set_option("--genotype", "%(IN_VCFFILE)s")
        params.set_option("--intervals", "%(IN_INTERVALS)s")

        params.set_kwargs(IN_VCFFILE=infile,
                          IN_TABIX=infile + ".tbi",
                          IN_INTERVALS=bedfile,
                          OUT_STDOUT=outfile)

        return {"command": params}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()
        description = "<BuildRegions: '%s' -> '%s'>" % (parameters.infile,
                                                        parameters.outfile)
        CommandNode.__init__(self,
                             description=description,
                             command=command,
                             dependencies=parameters.dependencies)


class SampleRegionsNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, infile, bedfile, outfile, dependencies=()):
        params = factory.new("sample_pileup")
        params.set_option("--genotype", "%(IN_PILEUP)s")
        params.set_option("--intervals", "%(IN_INTERVALS)s")
        params.set_kwargs(IN_PILEUP=infile,
                          IN_INTERVALS=bedfile,
                          OUT_STDOUT=outfile)

        return {"command": params}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()

        description = "<SampleRegions: '%s' -> '%s'>" \
            % (parameters.infile, parameters.outfile)

        CommandNode.__init__(self,
                             description=description,
                             command=command,
                             dependencies=parameters.dependencies)
