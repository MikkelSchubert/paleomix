#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
"""
SAMTools - Tools (written in C using htslib) for manipulating next-generation
sequencing data

https://github.com/samtools/samtools
"""

from __future__ import annotations

import os
from typing import Iterable, Literal

from paleomix.common import versions
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
        if options is None:
            options = {}

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
        index_format: Literal[".bai", ".csi"] = ".bai",
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        command = AtomicCmd(
            ["samtools", "index"],
            requirements=[SAMTOOLS_VERSION],
        )

        if index_format == ".csi":
            command.append("-c")
        elif index_format != ".bai":
            raise ValueError(f"Unknown BAM index format {index_format!r}")

        command.append_options(options)
        command.append(
            InputFile(infile),
            OutputFile(infile + index_format),
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"creating {index_format[1:].upper()} index for {infile}",
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
        index_format: str = ".bai",
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        if method not in ("stats", "idxstats", "flagstats"):
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
        index_format: Literal[".bai", ".csi"] | None = None,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        in_files = tuple(in_files)
        if not in_files:
            raise ValueError("no input files for samtools merge")
        elif len(in_files) == 1:
            # FIXME: hardlinking is faster, but could have unintended side-effects
            cmd = AtomicCmd(["cp", InputFile(in_files[0]), OutputFile(out_file)])
            self.index = None

            threads = 1
        else:
            cmd = AtomicCmd(["samtools", "merge"])
            cmd.append_options(options)

            self.index, requirement = _try_write_index(out_file, index_format)
            if self.index is not None:
                cmd.add_extra_files([OutputFile(self.index), OutputFile(out_file)])
                cmd.append("--write-index")
                # The location of the output index (and its format) can be specified by
                # writing the output file as "${outfile}##idx##${outindex}"
                cmd.append(
                    "%(TEMP_DIR)s/{}##idx##%(TEMP_DIR)s/{}".format(
                        os.path.basename(out_file),
                        os.path.basename(self.index),
                    )
                )
            else:
                cmd.append(OutputFile(out_file))

            for in_file in in_files:
                cmd.append(InputFile(in_file))

            cmd.add_requirement(requirement)

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
        *,
        in_bams: Iterable[str],
        out_bam: str,
        out_stats: str | None = None,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

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


def _get_number_of_threads(options: OptionsType, default: int = 1) -> int:
    if "-@" in options and "--threads" in options:
        raise ValueError(f"cannot use both -@ and --threads: {options!r}")

    value = options.get("-@", options.get("--threads", default))
    if not isinstance(value, int):
        raise TypeError(f"invalid number of samtools threads: {value}")

    # -@/--threads specify threads in *addition* to the main thread, but in practice
    # the number of cores used seems to be closer to the that value and not value + 1
    return max(1, value)


def _try_write_index(
    out_bam: str,
    index_format: str | None,
) -> tuple[str | None, versions.Requirement]:
    if index_format not in (None, ".bai", ".csi"):
        raise ValueError(index_format)

    try:
        if index_format is not None and SAMTOOLS_VERSION_1_10.check():
            return out_bam + index_format, SAMTOOLS_VERSION_1_10
    except versions.RequirementError:
        pass

    return None, SAMTOOLS_VERSION
