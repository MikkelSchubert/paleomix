#!/usr/bin/env python3
"""
SAMTools - Tools (written in C using htslib) for manipulating next-generation
sequencing data

https://github.com/samtools/samtools
"""
import os
from typing import Iterable, Optional, Tuple

import paleomix.common.versions as versions
from paleomix.common.command import (
    AtomicCmd,
    InputFile,
    OptionsType,
    OutputFile,
    ParallelCmds,
    SequentialCmds,
    TempInputFile,
    TempOutputFile,
)
from paleomix.common.fileutils import describe_files
from paleomix.node import CommandNode, Node

_VERSION_REGEX = r"Version: (\d+\.\d+)(?:\.(\d+))?"

SAMTOOLS_VERSION = versions.Requirement(
    call=("samtools",),
    regexp=_VERSION_REGEX,
    specifiers=">=1.6.0",
)

# Version required for --write-index
SAMTOOLS_VERSION_1_10 = versions.Requirement(
    call=("samtools",),
    regexp=_VERSION_REGEX,
    specifiers=">=1.10.0",
)

BCFTOOLS_VERSION = versions.Requirement(
    call=("bcftools",),
    regexp=_VERSION_REGEX,
    specifiers=">=1.4.0",
)

TABIX_VERSION = versions.Requirement(
    call=("tabix",),
    regexp=_VERSION_REGEX,
    specifiers=">=1.3.1",
)


class TabixIndexNode(CommandNode):
    """Tabix indexes a BGZip compressed VCF or pileup file."""

    def __init__(
        self,
        infile: str,
        preset: str = "vcf",
        options: OptionsType = {},
        dependencies: Iterable[Node] = (),
    ):
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

    def __init__(self, infile: str, dependencies: Iterable[Node] = ()):
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

    def __init__(
        self,
        infile: str,
        index_format: str = ".bai",
        options: OptionsType = {},
        dependencies: Iterable[Node] = (),
    ):
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
        method: str,
        infile: str,
        outfile: str,
        index_format: str = ".bai",
        options: OptionsType = {},
        dependencies: Iterable[Node] = (),
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
    def __init__(
        self,
        in_files: Iterable[str],
        out_file: str,
        index_format: Optional[str] = None,
        options: OptionsType = {},
        dependencies: Iterable[Node] = (),
    ):
        in_files = tuple(in_files)
        if not in_files:
            raise ValueError("no input files for samtools merge")
        elif len(in_files) == 1:
            # FIXME: hardlinking is faster, but could have unintended side-effects
            cmd = AtomicCmd(["cp", InputFile(in_files[0]), OutputFile(out_file)])
            self.index = None

            threads = 1
        else:
            self.index, requirement = _try_write_index(out_file, index_format)

            cmd = AtomicCmd(
                ["samtools", "merge"],
                requirements=[requirement],
            )

            if self.index is not None:
                cmd.append("--write-index")
                cmd.add_extra_files([OutputFile(self.index)])

            cmd.append_options(options)
            cmd.append(OutputFile(out_file))
            for in_file in in_files:
                cmd.append(InputFile(in_file))

            threads = _get_number_of_threads(options)

        CommandNode.__init__(
            self,
            command=cmd,
            description="merging %i files into %s" % (len(in_files), out_file),
            threads=threads,
            dependencies=dependencies,
        )


class MarkDupNode(CommandNode):
    def __init__(
        self,
        in_bams: Iterable[str],
        out_bam: str,
        out_stats: Optional[str] = None,
        options: OptionsType = {},
        dependencies: Iterable[Node] = (),
    ):
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

        fixed_options: OptionsType = {"-T": "%(TEMP_DIR)s/markdup"}
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


def merge_bam_files_command(input_files: Iterable[str]):
    merge = AtomicCmd(
        ["samtools", "merge", "-u", "-"],
        stdout=AtomicCmd.PIPE,
        requirements=[SAMTOOLS_VERSION],
    )

    for filename in input_files:
        merge.append(InputFile(filename))

    return merge


def _get_number_of_threads(options: OptionsType, default: int = 1) -> int:
    if "-@" in options and "--threads" in options:
        raise ValueError("cannot use both -@ and --threads: {!r}".format(options))

    value = options.get("-@", options.get("--threads", default))
    if not isinstance(value, int):
        raise ValueError(f"invalid number of samtools threads: {value}")

    # -@/--threads specify threads in *addition* to the main thread, but in practice
    # the number of cores used seems to be closer to the that value and not value + 1
    return max(1, value)


def _try_write_index(
    out_bam: str,
    index_format: Optional[str],
) -> Tuple[Optional[str], versions.Requirement]:
    if index_format not in (None, ".bai", ".csi"):
        raise ValueError(index_format)

    try:
        if index_format is not None and SAMTOOLS_VERSION_1_10.check():
            return out_bam + index_format, SAMTOOLS_VERSION_1_10
    except versions.RequirementError:
        pass

    return None, SAMTOOLS_VERSION
