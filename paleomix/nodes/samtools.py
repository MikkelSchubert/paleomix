#!/usr/bin/env python3
"""
SAMTools - Tools (written in C using htslib) for manipulating next-generation
sequencing data

https://github.com/samtools/samtools
"""
import os

from paleomix.node import CommandNode
from paleomix.atomiccmd.builder import AtomicCmdBuilder
from paleomix.atomiccmd.command import AtomicCmd

from paleomix.common.fileutils import reroot_path
import paleomix.common.versions as versions


_VERSION_REGEX = r"Version: (\d+)\.(\d+)(?:\.(\d+))?"

SAMTOOLS_VERSION = versions.Requirement(
    call=("samtools",), search=_VERSION_REGEX, checks=versions.GE(1, 3, 1)
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
            assert False, "Unxpected preset: %r" % preset

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


class BAMStatsNode(CommandNode):
    METHODS = ("stats", "idxstats", "flagstats")

    def __init__(self, method, infile, outfile, index_format=".bai", dependencies=()):
        if method not in self.METHODS:
            raise ValueError(method)

        command = AtomicCmd(
            ["samtools", method, "%(IN_BAM)s"],
            IN_BAM=infile,
            IN_IDX=infile + index_format if method == "idxstats" else None,
            OUT_STDOUT=outfile,
            CHECK_SAM=SAMTOOLS_VERSION,
        )

        CommandNode.__init__(
            self,
            description="collecting %s for %s" % (method, infile),
            command=command,
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
