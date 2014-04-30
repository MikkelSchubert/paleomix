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
    $ palemix [...]
"""
from pypeline.node import \
    CommandNode, \
    Node
from pypeline.atomiccmd.command import \
    AtomicCmd
from pypeline.atomiccmd.sets import \
    ParallelCmds
from pypeline.common.fileutils import \
    describe_files
from pypeline.nodes.picard import \
    MultiBAMInput, \
    MultiBAMInputNode
from pypeline.atomiccmd.builder import \
    AtomicCmdBuilder, \
    create_customizable_cli_parameters, \
    use_customizable_cli_parameters
from pypeline.common.fileutils import \
    reroot_path, \
    move_file

import pypeline.tools.bam_stats.coverage as coverage
import pypeline.tools.factory as factory


class DuplicateHistogramNode(MultiBAMInputNode):
    """Node for calling the 'paleomix duphist' command.

    Takes 1 or more BAMs as imput, requiring a config object in order to run
    MergeSamFiles.jar to merge these files. The output is a histogram of PCR
    duplicate counts, usable as input for the 'preseq' tool.
    """

    def __init__(self, config, input_files, output_file, dependencies=()):
        bam_input = MultiBAMInput(config, input_files)
        duphist_command = factory.new("duphist")
        duphist_command.add_value('%(TEMP_IN_BAM)s')
        duphist_command.set_kwargs(TEMP_IN_BAM=bam_input.pipe,
                                   OUT_STDOUT=output_file)
        duphist_command = duphist_command.finalize()

        commands = ParallelCmds(bam_input.commands + [duphist_command])

        description = "<DuplicateHistogram: %s -> %r>" \
            % (describe_files(input_files), output_file)
        MultiBAMInputNode.__init__(self,
                                   bam_input=bam_input,
                                   command=commands,
                                   description=description,
                                   dependencies=dependencies)


class CoverageNode(MultiBAMInputNode):
    def __init__(self, config, target_name, input_files, output_file,
                 regions_file=None, dependencies=()):
        bam_input = MultiBAMInput(config, input_files)
        if len(bam_input.files) > 1 and regions_file:
            raise ValueError("Coverage for regions require single, "
                             "indexed input BAM file.")

        builder = factory.new("coverage")
        builder.add_value("%(TEMP_IN_BAM)s")
        builder.add_value("%(OUT_FILE)s")
        builder.set_option("--target-name", target_name)
        builder.set_kwargs(TEMP_IN_BAM=bam_input.pipe,
                           OUT_FILE=output_file)

        if regions_file:
            builder.set_option('--regions-file', '%(IN_REGIONS)s')
            builder.set_kwargs(IN_REGIONS=regions_file)

        command = ParallelCmds(bam_input.commands + [builder.finalize()])
        description = "<Coverage: %s -> '%s'>" \
            % (describe_files(bam_input.files), output_file)
        MultiBAMInputNode.__init__(self,
                                   bam_input=bam_input,
                                   command=command,
                                   description=description,
                                   dependencies=dependencies)


class MergeCoverageNode(Node):
    def __init__(self, input_files, output_file, dependencies=()):
        self._output_file = output_file

        Node.__init__(self,
                      description="<MergeCoverage: '%s' -> '%s'>"
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
        bam_input = MultiBAMInput(config, input_files)
        if len(bam_input.files) > 1 and regions_file:
            raise ValueError("DepthHistogram for regions require single, "
                             "indexed input BAM file.")

        builder = factory.new("depths")
        builder.add_value("%(TEMP_IN_BAM)s")
        builder.add_value("%(OUT_FILE)s")
        builder.set_option("--target-name", target_name)
        builder.set_kwargs(TEMP_IN_BAM=bam_input.pipe,
                           OUT_FILE=output_file)

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
        bam_input = MultiBAMInput(config, input_bams)

        builder = factory.new("rmdup_collapsed")
        builder.add_value("%(TEMP_IN_BAM)s")
        builder.set_kwargs(TEMP_IN_BAM=bam_input.pipe,
                           OUT_STDOUT=output_bam)
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
    @create_customizable_cli_parameters
    def customize(cls, reference, in_bam, in_vcf, outfile, dependencies=()):
        cat = factory.new("cat")
        cat.add_value("%(IN_VCF)s")
        cat.set_kwargs(IN_VCF=in_vcf,
                       OUT_STDOUT=AtomicCmd.PIPE)

        pileup = factory.new("create_pileup", outfile="%(OUT_PILEUP)s")
        pileup.add_value("%(IN_BAM)s")
        pileup.set_option("-f", "%(IN_REF)s")
        pileup.set_kwargs(IN_REF=reference,
                          IN_BAM=in_bam,
                          IN_STDIN=cat,
                          OUT_PILEUP=outfile,
                          OUT_TBI=outfile + ".tbi")

        return {"commands": {"cat": cat,
                             "pileup": pileup}}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        commands = parameters.commands
        commands = [commands[key].finalize() for key in ("cat", "pileup")]
        description = "<VCFPileup: '%s' -> '%s'>" % (parameters.in_bam,
                                                     parameters.outfile)
        CommandNode.__init__(self,
                             description=description,
                             command=ParallelCmds(commands),
                             dependencies=parameters.dependencies)


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