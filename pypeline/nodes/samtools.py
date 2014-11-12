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

from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd

from pypeline.common.fileutils import reroot_path, swap_ext
import pypeline.common.versions as versions


_VERSION_REGEX = r"Version: (\d+)\.(\d+)(?:\.(\d+))?"

# v0.2.0 was the pre-release version of v1.0, and lacks required features
_COMMON_CHECK = versions.And(versions.GE(0, 1, 18),
                             versions.LT(0, 2, 0))

SAMTOOLS_VERSION = versions.Requirement(call=("samtools",),
                                        search=_VERSION_REGEX,
                                        checks=_COMMON_CHECK)

BCFTOOLS_VERSION \
    = versions.Requirement(call=("bcftools",),
                           search=_VERSION_REGEX,
                           checks=_COMMON_CHECK)

TABIX_VERSION = versions.Requirement(call=("tabix",),
                                     search=_VERSION_REGEX,
                                     checks=versions.GE(0, 2, 5))


def samtools_compatible_wbu_mode():
    """Returns a writing mode for Pysam compatible with the current version of
    samtools; uncompressed output from Pysam 0.8.x cannot be read by older
    versions of samtools:

    https://github.com/pysam-developers/pysam/issues/43
    """
    import pysam

    if map(int, pysam.__version__.split(".")) < [0, 8, 0]:
        return "wbu"
    else:
        return "wb"


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
        cmd_index = AtomicCmd(["samtools", "index", "%(IN_BAM)s",
                               "%(OUT_BAI)s"],
                              IN_BAM=infile,
                              OUT_BAI=swap_ext(infile, ".bai"),
                              CHECK_SAM=SAMTOOLS_VERSION)

        CommandNode.__init__(self,
                             description="<BAMIndex: '%s'>" % (infile,),
                             command=cmd_index,
                             dependencies=dependencies)
