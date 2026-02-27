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

import paleomix.common.versions as versions
from paleomix.atomiccmd.builder import AtomicCmdBuilder
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import ParallelCmds
from paleomix.common.fileutils import describe_files, reroot_path
from paleomix.node import CommandNode

_VERSION_REGEX = r"Version: (\d+)\.(\d+)(?:\.(\d+))?"

SAMTOOLS_VERSION = versions.Requirement(
    call=("samtools",), search=_VERSION_REGEX, checks=versions.GE(1, 6, 0)
)

BCFTOOLS_VERSION = versions.Requirement(
    call=("bcftools",), search=_VERSION_REGEX, checks=versions.GE(1, 4, 0)
)

TABIX_VERSION = versions.Requirement(
    call=("tabix",), search=_VERSION_REGEX, checks=versions.GE(1, 3, 1)
)


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
            assert False, "Unexpected preset: %r" % preset

        self._infile = infile
        cmd_tabix = AtomicCmd(
            call + ["%(TEMP_IN_VCFFILE)s"],
            TEMP_IN_VCFFILE=os.path.basename(infile),
            IN_VCFFILE=infile,
            OUT_TBI=infile + ".tbi",
            CHECK_TABIX=TABIX_VERSION,
        )

        CommandNode.__init__(
            self,
            description="creating tabix %s index for %s" % (preset, infile),
            command=cmd_tabix,
            dependencies=dependencies,
        )

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
        cmd_faidx = AtomicCmd(
            ["samtools", "faidx", "%(TEMP_IN_FASTA)s"],
            TEMP_IN_FASTA=os.path.basename(infile),
            IN_FASTA=infile,
            OUT_TBI=infile + ".fai",
            CHECK_SAM=SAMTOOLS_VERSION,
        )

        CommandNode.__init__(
            self,
            description="creating FAI index for %s" % (infile,),
            command=cmd_faidx,
            dependencies=dependencies,
        )

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

    def __init__(self, infile, index_format=".bai", dependencies=()):
        if index_format == ".bai":
            samtools_call = ["samtools", "index", "%(IN_BAM)s", "%(OUT_IDX)s"]
        elif index_format == ".csi":
            samtools_call = ["samtools", "index", "-c", "%(IN_BAM)s", "%(OUT_IDX)s"]
        else:
            raise ValueError(
                "Unknown format type %r; expected .bai or .csi" % (index_format,)
            )

        command = AtomicCmd(
            samtools_call,
            IN_BAM=infile,
            OUT_IDX=infile + index_format,
            CHECK_SAM=SAMTOOLS_VERSION,
        )

        CommandNode.__init__(
            self,
            description="creating %s index for %s" % (index_format[1:].upper(), infile),
            command=command,
            dependencies=dependencies,
        )


class MarkDupNode(CommandNode):
    def __init__(
        self,
        in_bams,
        out_bam,
        out_stats=None,
        keep_duplicates=True,
        dependencies=(),
    ):
        in_bams = tuple(in_bams)
        if len(in_bams) > 1:
            merge = merge_bam_files_command(in_bams)

            markdup = AtomicCmdBuilder(
                ["samtools", "markdup", "-", "%(OUT_BAM)s"],
                IN_STDIN=merge,
                OUT_BAM=out_bam,
                # Stderr is piped instead of saved using -f to support samtools < v1.10
                OUT_STDERR=out_stats,
                CHECK_VERSION=SAMTOOLS_VERSION,
                set_cwd=True,
            )

            threads = 2
        elif in_bams:
            merge = None
            markdup = AtomicCmdBuilder(
                ["samtools", "markdup", "%(IN_BAM)s", "%(OUT_BAM)s"],
                IN_BAM=in_bams[0],
                OUT_BAM=out_bam,
                OUT_STDERR=out_stats,
                CHECK_VERSION=SAMTOOLS_VERSION,
                set_cwd=True,
            )

            threads = 1
        else:
            raise ValueError(in_bams)

        if out_stats is not None:
            markdup.add_option("-s", None)

        if not keep_duplicates:
            markdup.add_option("-r", None)

        markdup = markdup.finalize()

        CommandNode.__init__(
            self,
            command=markdup if merge is None else ParallelCmds([merge, markdup]),
            description="marking PCR duplicates in {}".format(describe_files(in_bams)),
            dependencies=dependencies,
            threads=threads,
        )


class BAMMergeNode(CommandNode):
    def __init__(self, input_bams, output_bam, dependencies=()):
        in_files = tuple(input_bams)
        if not in_files:
            raise ValueError("no input files for samtools merge")
        elif len(in_files) == 1:
            # FIXME: hardlinking is faster, but could have unintended side-effects
            cmd = AtomicCmd(
                ["cp", "%(IN_BAM)s", "%(OUT_BAM)s"],
                IN_BAM=in_files[0],
                OUT_BAM=output_bam,
            )
        else:
            cmd = AtomicCmdBuilder(
                ["samtools", "merge", "%(OUT_BAM)s"],
                OUT_BAM=output_bam,
                CHECK_VERSION=SAMTOOLS_VERSION,
            )

            cmd.add_multiple_values(in_files)
            cmd = cmd.finalize()

        CommandNode.__init__(
            self,
            command=cmd,
            description=f"merging {len(in_files)} files into {output_bam}",
            dependencies=dependencies,
        )


def merge_bam_files_command(input_files):
    merge = AtomicCmdBuilder(
        ["samtools", "merge", "-u", "-"],
        OUT_STDOUT=AtomicCmd.PIPE,
        CHECK_VERSION=SAMTOOLS_VERSION,
    )

    merge.add_multiple_values(input_files)

    return merge.finalize()
