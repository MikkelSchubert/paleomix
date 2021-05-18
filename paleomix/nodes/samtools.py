#!/usr/bin/env python3
"""
SAMTools - Tools (written in C using htslib) for manipulating next-generation
sequencing data

https://github.com/samtools/samtools
"""
import os
import warnings

from paleomix.node import CommandNode
from paleomix.atomiccmd.command import (
    AtomicCmd,
    InputFile,
    OutputFile,
    TempInputFile,
    TempOutputFile,
)
from paleomix.atomiccmd.sets import ParallelCmds, SequentialCmds
from paleomix.common.fileutils import describe_files

import paleomix.common.versions as versions


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
    """Tabix indexes a BGZip compressed VCF or pileup file."""

    def __init__(self, infile, preset="vcf", options={}, dependencies=()):
        if preset not in ("vcf", "gff", "bed", "sam"):
            raise ValueError(preset)

        basename = os.path.basename(infile)
        infile = os.path.abspath(infile)

        # Tabix does not support a custom output path, so we create a symlink to the
        # input file in the temporary folder and index that.
        link = AtomicCmd(["ln", "-s", InputFile(infile), TempOutputFile(basename)])

        tabix = AtomicCmd(
            ["tabix", "-p", preset, TempInputFile(basename)],
            extra_files=[OutputFile(infile + ".tbi")],
            requirements=[TABIX_VERSION],
        )

        tabix.append_options(options)

        CommandNode.__init__(
            self,
            command=SequentialCmds([link, tabix]),
            description="creating tabix %s index for %s" % (preset, infile),
            dependencies=dependencies,
        )


class FastaIndexNode(CommandNode):
    """Indexed a FASTA file using 'samtools faidx'."""

    def __init__(self, infile, dependencies=()):
        basename = os.path.basename(infile)

        # faidx does not support a custom output path, so we create a symlink to the
        # input file in the temporary folder and index that.
        link = AtomicCmd(
            [
                "ln",
                "-s",
                InputFile(os.path.abspath(infile)),
                TempOutputFile(basename),
            ]
        )

        faidx = AtomicCmd(
            ["samtools", "faidx", TempInputFile(basename)],
            extra_files=[OutputFile(infile + ".fai")],
            requirements=[SAMTOOLS_VERSION],
        )

        CommandNode.__init__(
            self,
            description="creating FAI index for %s" % (infile,),
            command=SequentialCmds([link, faidx]),
            dependencies=dependencies,
        )


class BAMIndexNode(CommandNode):
    """Indexed a BAM file using 'samtools index'."""

    def __init__(self, infile, index_format=".bai", options={}, dependencies=()):
        command = AtomicCmd(
            ["samtools", "index"],
            requirements=[SAMTOOLS_VERSION],
        )

        if index_format == ".csi":
            command.append("-c")
        elif index_format != ".bai":
            raise ValueError("Unknown BAM index format %r" % (index_format,))

        command.append_options(options)
        command.append(
            InputFile(infile),
            OutputFile(infile + index_format),
        )

        CommandNode.__init__(
            self,
            command=command,
            description="creating %s index for %s" % (index_format[1:].upper(), infile),
            threads=_get_number_of_threads(options),
            dependencies=dependencies,
        )


class BAMStatsNode(CommandNode):
    METHODS = ("stats", "idxstats", "flagstats")

    def __init__(
        self,
        method,
        infile,
        outfile,
        index_format=".bai",
        options={},
        dependencies=(),
    ):
        if method not in self.METHODS:
            raise ValueError(method)

        command = AtomicCmd(
            ["samtools", method, InputFile(infile)],
            stdout=outfile,
            requirements=[SAMTOOLS_VERSION],
        )

        command.append_options(options)

        if method == "idxstats":
            command.add_extra_files([InputFile(infile + index_format)])

        CommandNode.__init__(
            self,
            command=command,
            description="collecting %s for %s" % (method, infile),
            threads=_get_number_of_threads(options),
            dependencies=dependencies,
        )


class BAMMergeNode(CommandNode):
    def __init__(self, in_files, out_file, options={}, dependencies=()):
        in_files = tuple(in_files)
        if len(in_files) <= 1:
            warnings.warn("creating {!r} from single input file".format(out_file))

        cmd = AtomicCmd(
            ["samtools", "merge"],
            requirements=[SAMTOOLS_VERSION],
        )

        cmd.append_options(options)
        cmd.append(OutputFile(out_file))
        for in_file in in_files:
            cmd.append(InputFile(in_file))

        CommandNode.__init__(
            self,
            command=cmd,
            description="merging %i files into %s" % (len(in_files), out_file),
            threads=_get_number_of_threads(options),
            dependencies=dependencies,
        )


class MarkDupNode(CommandNode):
    def __init__(self, in_bams, out_bam, out_stats=None, options={}, dependencies=()):
        in_bams = tuple(in_bams)
        if len(in_bams) > 1:
            merge = AtomicCmd(
                ["samtools", "merge", "-u", "-"],
                stdout=AtomicCmd.PIPE,
                requirements=[SAMTOOLS_VERSION],
            )

            for in_file in in_bams:
                merge.append(InputFile(in_file))

            markdup = AtomicCmd(
                ["samtools", "markdup", "-", OutputFile(out_bam)],
                stdin=merge,
                # Stderr is piped instead of saved using -f to support samtools < v1.10
                stderr=out_stats,
                requirements=[SAMTOOLS_VERSION],
            )

            command = ParallelCmds([merge, markdup])
        else:
            (in_file,) = in_bams

            command = markdup = AtomicCmd(
                ["samtools", "markdup", InputFile(in_file), OutputFile(out_bam)],
                requirements=[SAMTOOLS_VERSION],
            )

        fixed_options = {"-T": "%(TEMP_DIR)s/markdup"}
        if out_stats is not None:
            fixed_options["-s"] = None
        markdup.merge_options(
            user_options=options,
            fixed_options=fixed_options,
            blacklisted_options=["-f"],
        )

        CommandNode.__init__(
            self,
            command=command,
            description="marking PCR duplicates in {}".format(describe_files(in_bams)),
            threads=_get_number_of_threads(options),
            dependencies=dependencies,
        )


def merge_bam_files_command(input_files):
    merge = AtomicCmd(
        ["samtools", "merge", "-u", "-"],
        stdout=AtomicCmd.PIPE,
        requirements=[SAMTOOLS_VERSION],
    )

    for filename in input_files:
        merge.append(InputFile(filename))

    return merge


def _get_number_of_threads(options, default=1):
    if "-@" in options and "--threads" in options:
        raise ValueError("cannot use both -@ and --threads: {!r}".format(options))

    # -@/--threads specify threads in *addition* to the main thread, but in practice
    # the number of cores used seems to be closer to the that value and not value + 1
    return max(1, options.get("-@", options.get("--threads", default)))
