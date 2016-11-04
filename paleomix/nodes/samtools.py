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
import os

from paleomix.node import CommandNode
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import SequentialCmds

from paleomix.common.fileutils import reroot_path, swap_ext
import paleomix.common.versions as versions


_VERSION_REGEX = r"Version: (\d+)\.(\d+)(?:\.(\d+))?"

# v0.2.0 was the pre-release version of v1.0, and lacks required features
_COMMON_CHECK = versions.Or(versions.EQ(0, 1, 19),
                            versions.GE(1, 0, 0))

SAMTOOLS_VERSION = versions.Requirement(call=("samtools",),
                                        search=_VERSION_REGEX,
                                        checks=_COMMON_CHECK)

SAMTOOLS_VERSION_0119 = versions.Requirement(call=("samtools",),
                                             search=_VERSION_REGEX,
                                             checks=versions.EQ(0, 1, 19))

BCFTOOLS_VERSION_0119 \
    = versions.Requirement(call=("bcftools",),
                           search=_VERSION_REGEX,
                           checks=versions.EQ(0, 1, 19))

TABIX_VERSION = versions.Requirement(call=("tabix",),
                                     search=_VERSION_REGEX,
                                     checks=versions.GE(0, 2, 5))


class TabixIndexNode(CommandNode):
    """Tabix indexes a BGZip compressed VCF or pileup file.

    The class currently supports the following presets:
        - vcf -- BGZipped VCF file.
        - pileup -- BGZipped pileup (non-binary) as produced by 'mpileup'.
    """

    def __init__(self, infile, preset="vcf", dependencies=()):
        assert infile.lower().endswith(".bgz")
        if preset == "pileup":
            call = ["tabix", "-s", 1, "-b", 2, "-e", 2]
        elif preset == "vcf":
            call = ["tabix", "-p", preset]
        else:
            assert False, "Unxpected preset: %r" % preset

        self._infile = infile
        cmd_tabix = AtomicCmd(call + ["%(TEMP_IN_VCFFILE)s"],
                              TEMP_IN_VCFFILE=os.path.basename(infile),
                              IN_VCFFILE=infile,
                              OUT_TBI=infile + ".tbi",
                              CHECK_TABIX=TABIX_VERSION)

        CommandNode.__init__(self,
                             description="<TabixIndex (%s): '%s'>" % (preset,
                                                                      infile,),
                             command=cmd_tabix,
                             dependencies=dependencies)

    def _setup(self, config, temp):
        """See CommandNode._setup."""
        infile = os.path.abspath(self._infile)
        outfile = reroot_path(temp, self._infile)
        os.symlink(infile, outfile)

        CommandNode._setup(self, config, temp)

    def _teardown(self, config, temp):
        """See CommandNode._teardown."""
        os.remove(reroot_path(temp, self._infile))

        CommandNode._teardown(self, config, temp)


class FastaIndexNode(CommandNode):
    """Indexed a FASTA file using 'samtools faidx'."""

    def __init__(self, infile, dependencies=()):
        self._infile = infile
        cmd_faidx = AtomicCmd(["samtools", "faidx", "%(TEMP_IN_FASTA)s"],
                              TEMP_IN_FASTA=os.path.basename(infile),
                              IN_FASTA=infile,
                              OUT_TBI=infile + ".fai",
                              CHECK_SAM=SAMTOOLS_VERSION)

        CommandNode.__init__(self,
                             description="<FastaIndex: '%s'>" % (infile,),
                             command=cmd_faidx,
                             dependencies=dependencies)

    def _setup(self, config, temp):
        """See CommandNode._setup."""
        infile = os.path.abspath(self._infile)
        outfile = reroot_path(temp, self._infile)
        os.symlink(infile, outfile)

        CommandNode._setup(self, config, temp)

    def _teardown(self, config, temp):
        """See CommandNode._teardown."""
        os.remove(reroot_path(temp, self._infile))

        CommandNode._teardown(self, config, temp)


class BAMIndexNode(CommandNode):
    """Indexed a BAM file using 'samtools index'."""

    def __init__(self, infile, dependencies=()):
        basename = os.path.basename(infile)

        cmd_link = AtomicCmd(["ln", "-s", "%(IN_BAM)s", "%(TEMP_OUT_BAM)s"],
                             IN_BAM=infile,
                             TEMP_OUT_BAM=basename,
                             set_cwd=True)

        cmd_index = AtomicCmd(["samtools", "index", "%(TEMP_IN_BAM)s"],
                              TEMP_IN_BAM=basename,
                              CHECK_SAM=SAMTOOLS_VERSION)

        cmd_rename = AtomicCmd(["mv", "%(TEMP_IN_BAM)s", "%(OUT_BAM)s"],
                               TEMP_IN_BAM=basename + ".bai",
                               OUT_BAM=swap_ext(infile, ".bai"))

        commands = SequentialCmds((cmd_link, cmd_index, cmd_rename))

        CommandNode.__init__(self,
                             description="<BAMIndex: '%s'>" % (infile,),
                             command=commands,
                             dependencies=dependencies)


class RMDuplicatesNode(CommandNode):
    """Remove PCR duplicates from BAM file."""

    def __init__(self, input_bam, output_bam, se_reads=False, force_se=False,
                 dependencies=()):
        call = ["samtools", "rmdup"]
        if se_reads:
            call.append("-s")
        if force_se:
            call.append("-S")

        command = AtomicCmd(call + ["%(IN_BAM)s", "%(OUT_BAM)s"],
                            IN_BAM=input_bam,
                            OUT_BAM=output_bam,
                            CHECK_SAM=SAMTOOLS_VERSION)

        CommandNode.__init__(self,
                             description="<Samtools rmdup: %r -> %r>"
                             % (input_bam, output_bam),
                             command=command,
                             dependencies=dependencies)


class FilterBAMNode(CommandNode):
    """Filter BAM file using samtools view."""

    def __init__(self, input_bam, output_bam, require_flags=0, exclude_flags=0,
                 dependencies=()):
        call = ["samtools", "view", "-b"]
        if require_flags:
            call.extend(("-f", hex(require_flags)))
        if exclude_flags:
            call.extend(("-F", hex(exclude_flags)))

        command = AtomicCmd(call + ["%(IN_BAM)s"],
                            IN_BAM=input_bam,
                            OUT_STDOUT=output_bam,
                            CHECK_SAM=SAMTOOLS_VERSION)

        CommandNode.__init__(self,
                             description="<SAMTools view: %r -> %r>"
                             % (input_bam, output_bam),
                             command=command,
                             dependencies=dependencies)
