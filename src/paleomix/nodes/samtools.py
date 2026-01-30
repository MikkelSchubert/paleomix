# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
"""
SAMTools - Tools (written in C using htslib) for manipulating next-generation
sequencing data

https://github.com/samtools/samtools
"""

from __future__ import annotations

import os
from collections.abc import Iterable
from typing import Literal

from paleomix.common import versions
from paleomix.common.bamfiles import get_idx_filename
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
    specifiers=">=1.10.0",
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
        *,
        infile: str,
        preset: Literal["vcf", "gff", "bed", "sam"] = "vcf",
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
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
            description=f"creating tabix {preset} index for {infile}",
            dependencies=dependencies,
        )


class FastaIndexNode(CommandNode):
    """Indexed a FASTA file using 'samtools faidx'."""

    def __init__(self, infile: str, dependencies: Iterable[Node] = ()) -> None:
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
            description=f"creating FAI index for {infile}",
            command=SequentialCmds([link, faidx]),
            dependencies=dependencies,
        )


class BAMIndexNode(CommandNode):
    """Indexed a BAM file using 'samtools index'."""

    def __init__(
        self,
        *,
        infile: str,
        bam_index: Literal[".bai", ".csi"] = ".bai",
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = AtomicCmd(
            ["samtools", "index"],
            requirements=[SAMTOOLS_VERSION],
        )

        if bam_index == ".csi":
            command.append("-c")
        elif bam_index != ".bai":
            raise ValueError(f"Unknown BAM index format {bam_index!r}")

        command.append_options(options)
        command.append(
            InputFile(infile),
            OutputFile(get_idx_filename(infile, bam_index)),
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"creating index for {infile}",
            threads=_get_number_of_threads(options),
            dependencies=dependencies,
        )


class BAMStatsNode(CommandNode):
    def __init__(
        self,
        *,
        method: Literal["stats", "idxstats", "flagstats"],
        infile: str,
        outfile: str,
        reference: str,
        bam_index: Literal[".bai", ".csi"] = ".bai",
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = AtomicCmd(
            ["samtools", method, InputFile(infile)],
            stdout=outfile,
            requirements=[SAMTOOLS_VERSION],
        )

        command.append_options(options)

        if method == "stats":
            command.append("--reference", InputFile(reference))
            command.add_extra_files([InputFile(f"{reference}.fai")])
        elif method == "idxstats":
            command.add_extra_files([InputFile(get_idx_filename(infile, bam_index))])
        elif method != "flagstats":
            raise ValueError(method)

        CommandNode.__init__(
            self,
            command=command,
            description=f"collecting {method} for {infile}",
            threads=_get_number_of_threads(options),
            dependencies=dependencies,
        )


class BAMMergeNode(CommandNode):
    def __init__(
        self,
        *,
        in_files: Iterable[str],
        out_file: str,
        bam_index: Literal[".bai", ".csi"] = ".bai",
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        in_files = tuple(in_files)
        if not in_files:
            raise ValueError("no input files for samtools merge")
        elif len(in_files) == 1:
            # FIXME: hardlinking is faster, but could have unintended side-effects
            cmd = AtomicCmd(["cp", InputFile(in_files[0]), OutputFile(out_file)])
            self.index = None

            threads = 1
        else:
            self.index = get_idx_filename(out_file, bam_index)

            cmd = AtomicCmd(
                ["samtools", "merge", "--write-index"],
                set_cwd=True,
                requirements=[SAMTOOLS_VERSION],
                extra_files=[OutputFile(self.index), OutputFile(out_file)],
            )
            cmd.append_options(options)

            # The location of the output index (and its format) can be specified by
            # writing the output file as "${outfile}##idx##${outindex}"
            cmd.append(
                f"{os.path.basename(out_file)}##idx##{os.path.basename(self.index)}"
            )

            for in_file in in_files:
                cmd.append(InputFile(in_file))

            threads = _get_number_of_threads(options)

        CommandNode.__init__(
            self,
            command=cmd,
            description=f"merging {len(in_files)} files into {out_file}",
            threads=threads,
            dependencies=dependencies,
        )


class BAMToCRAMNode(CommandNode):
    def __init__(
        self,
        *,
        in_bam: str,
        out_cram: str,
        reference: str,
        threads: int,
        dependencies: Iterable[Node] = (),
    ) -> None:
        command = AtomicCmd(
            [
                "samtools",
                "view",
                "--cram",
                "--write-index",
                "--threads",
                threads,
                "--reference",
                InputFile(reference),
                "--output",
                OutputFile(out_cram),
                InputFile(in_bam),
            ],
            extra_files=[
                OutputFile(f"{out_cram}.crai"),
            ],
            requirements=[SAMTOOLS_VERSION],
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"converting {in_bam} to {out_cram}",
            threads=max(1, threads),
            dependencies=dependencies,
        )


class MarkDupNode(CommandNode):
    def __init__(
        self,
        *,
        in_bams: Iterable[str],
        out_bam: str,
        out_stats: str | None = None,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
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
            description=f"marking PCR duplicates in {describe_files(in_bams)}",
            threads=_get_number_of_threads(options),
            dependencies=dependencies,
        )


def merge_bam_files_command(input_files: Iterable[str]) -> AtomicCmd:
    merge = AtomicCmd(
        ["samtools", "merge", "-u", "-"],
        stdout=AtomicCmd.PIPE,
        requirements=[SAMTOOLS_VERSION],
    )

    for filename in input_files:
        merge.append(InputFile(filename))

    return merge


def _get_number_of_threads(options: OptionsType | None, default: int = 1) -> int:
    if options is None:
        return default
    elif "-@" in options and "--threads" in options:
        raise ValueError(f"cannot use both -@ and --threads: {options!r}")

    value = options.get("-@", options.get("--threads", default))
    if not isinstance(value, int):
        raise TypeError(f"invalid number of samtools threads: {value}")

    # -@/--threads specify threads in *addition* to the main thread, but in practice
    # the number of cores used seems to be closer to the that value and not value + 1
    return max(1, value)
