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
from __future__ import annotations

import functools
import os
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

from paleomix.common import versions
from paleomix.common.command import (
    AtomicCmd,
    AtomicFileTypes,
    Executable,
    InputFile,
    OptionsType,
    OutputFile,
    ParallelCmds,
    TempOutputFile,
)
from paleomix.common.fileutils import PathTypes, describe_paired_files
from paleomix.node import CommandNode, Node
from paleomix.nodes.samtools import SAMTOOLS_VERSION
from paleomix.tools import factory

# Index files used by BWA and BWA-MEM2 respectively
BWA_INDEX_EXT = (".amb", ".ann", ".bwt", ".pac", ".sa")
BWA_MEM2_INDEX_EXT = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")

BWA_VERSION = versions.Requirement(
    name="BWA",
    call=("bwa",),
    regexp=r"Version: (\d+\.\d+\.\d+)",
    specifiers=">=0.7.15",
)


BWA_MEM2_VERSION = versions.Requirement(
    name="BWA MEM2",
    call=("bwa-mem2", "version"),
    regexp=r"(\d+\.\d+\.\d+)",
    specifiers=">=2.2.1",
)


class BWAIndexNode(CommandNode):
    def __init__(self, *, input_file: str, dependencies: Iterable[Node] = ()) -> None:
        command = _new_bwa_command(
            (
                "bwa",
                "index",
                InputFile(input_file),
                "-p",
                TempOutputFile(input_file),
            ),
            reference=input_file,
            index_iotype=OutputFile,
        )

        CommandNode.__init__(
            self,
            command=command,
            description=f"creating BWA index for {input_file}",
            dependencies=dependencies,
        )


class BWAMem2IndexNode(CommandNode):
    def __init__(self, *, input_file: str, dependencies: Iterable[Node] = ()) -> None:
        in_file = InputFile(input_file)
        tmp_file = TempOutputFile(input_file)
        index_call = ("bwa-mem2", "index", in_file, "-p", tmp_file)

        CommandNode.__init__(
            self,
            command=_new_bwa_command(
                index_call,
                reference=input_file,
                index_ext=BWA_MEM2_INDEX_EXT,
                index_iotype=OutputFile,
                requirements=[BWA_MEM2_VERSION],
            ),
            description=f"creating BWA MEM2 index for {input_file}",
            dependencies=dependencies,
        )


class BWABacktrack(CommandNode):
    def __init__(
        self,
        *,
        input_file: str,
        output_file: str,
        reference: str,
        threads: int = 1,
        mapping_options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        threads = get_max_threads(reference, threads)
        command = _new_bwa_command(
            ("bwa", "aln", reference, InputFile(input_file)),
            reference=reference,
            stdout=output_file,
        )

        command.merge_options(
            user_options=mapping_options,
            fixed_options={"-t": threads},
        )

        CommandNode.__init__(
            self,
            command=command,
            description=get_node_description(
                name="BWA backtrack",
                input_files_1=input_file,
                reference=reference,
            ),
            threads=threads,
            dependencies=dependencies,
        )


class BWASamse(CommandNode):
    def __init__(
        self,
        *,
        input_file_fq: str,
        input_file_sai: str,
        output_file: str,
        reference: str,
        threads: int = 1,
        mapping_options: OptionsType | None = None,
        cleanup_options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if cleanup_options is None:
            cleanup_options = {}
        if mapping_options is None:
            mapping_options = {}

        samse = _new_bwa_command(
            (
                "bwa",
                "samse",
                reference,
                InputFile(input_file_sai),
                InputFile(input_file_fq),
            ),
            reference=reference,
            stdout=AtomicCmd.PIPE,
        )

        samse.append_options(mapping_options)

        cleanup = new_cleanup_command(
            stdin=samse,
            in_reference=reference,
            out_bam=output_file,
            max_threads=threads,
            options=cleanup_options,
        )

        CommandNode.__init__(
            self,
            command=ParallelCmds([samse, cleanup]),
            description=get_node_description(
                name="BWA samse",
                input_files_1=input_file_fq,
                reference=reference,
            ),
            threads=threads,
            dependencies=dependencies,
        )


class BWASampe(CommandNode):
    def __init__(
        self,
        *,
        input_file_fq_1: str,
        input_file_fq_2: str,
        input_file_sai_1: str,
        input_file_sai_2: str,
        output_file: str,
        reference: str,
        threads: int = 1,
        mapping_options: OptionsType | None = None,
        cleanup_options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if cleanup_options is None:
            cleanup_options = {}
        if mapping_options is None:
            mapping_options = {}

        sampe = _new_bwa_command(
            (
                "bwa",
                "sampe",
                reference,
                InputFile(input_file_sai_1),
                InputFile(input_file_sai_2),
                InputFile(input_file_fq_1),
                InputFile(input_file_fq_2),
            ),
            reference=reference,
            stdout=AtomicCmd.PIPE,
        )

        sampe.append_options(mapping_options)

        cleanup = new_cleanup_command(
            stdin=sampe,
            in_reference=reference,
            out_bam=output_file,
            max_threads=threads,
            paired_end=True,
            options=cleanup_options,
        )

        CommandNode.__init__(
            self,
            command=ParallelCmds([sampe, cleanup]),
            description=get_node_description(
                name="BWA sampe",
                input_files_1=input_file_fq_1,
                input_files_2=input_file_fq_2,
                reference=reference,
            ),
            threads=threads,
            dependencies=dependencies,
        )


class BWAAlgorithmNode(CommandNode):
    def __init__(
        self,
        *,
        input_file_1: str,
        output_file: str,
        reference: str,
        input_file_2: str | None = None,
        threads: int = 1,
        algorithm: Literal["mem", "bwasw", "mem2"] = "mem",
        alt_aware: bool = False,
        alt_optimize: bool = False,
        mapping_options: OptionsType | None = None,
        cleanup_options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if cleanup_options is None:
            cleanup_options = {}
        if mapping_options is None:
            mapping_options = {}

        if algorithm in ("mem", "bwasw"):
            aln_call = ("bwa", algorithm)
            index_ext = BWA_INDEX_EXT
            requirements = [BWA_VERSION]
        elif algorithm == "mem2":
            aln_call = ("bwa-mem2", "mem")
            index_ext = BWA_MEM2_INDEX_EXT
            requirements = [BWA_MEM2_VERSION]
        else:
            raise NotImplementedError(f"BWA algorithm {algorithm!r} not implemented")

        aln = _new_bwa_command(
            aln_call,
            reference=reference,
            index_ext=index_ext,
            requirements=requirements,
            stdout=AtomicCmd.PIPE,
        )

        if alt_aware or alt_optimize:
            if algorithm == "bwasw":
                raise NotImplementedError("bwasw is not ALT aware")

            aln.add_extra_files([InputFile(reference + ".alt")])

        threads = get_max_threads(reference, threads)

        aln.merge_options(
            user_options=mapping_options,
            fixed_options={
                "-t": threads,
                # Mark alternative hits as secondary; required by e.g. Picard
                "-M": None,
            },
        )

        # Positional arguments must be last on OSX
        aln.append(reference, InputFile(input_file_1))
        if input_file_2:
            aln.append(InputFile(input_file_2))

        cleanup = new_cleanup_command(
            stdin=aln,
            in_reference=reference,
            out_bam=output_file,
            max_threads=threads,
            paired_end=bool(input_file_1 and input_file_2),
            alt_aware=alt_aware,
            alt_optimize=alt_optimize,
            options=cleanup_options,
        )

        CommandNode.__init__(
            self,
            command=ParallelCmds([aln, cleanup]),
            description=get_node_description(
                name=f"BWA {algorithm}",
                input_files_1=input_file_1,
                input_files_2=input_file_2,
                reference=reference,
            ),
            threads=threads,
            dependencies=dependencies,
        )


def new_cleanup_command(
    *,
    stdin: AtomicCmd,
    in_reference: str,
    out_bam: str,
    max_threads: int = 1,
    paired_end: bool = False,
    alt_aware: bool = False,
    alt_optimize: bool = False,
    options: OptionsType | None = None,
) -> AtomicCmd:
    convert = factory.new(
        "cleanup",
        stdin=stdin,
        stdout=out_bam,
        requirements=[
            SAMTOOLS_VERSION,
        ],
    )

    fixed_options: OptionsType = {
        "--fasta": InputFile(in_reference),
        "--temp-prefix": TempOutputFile("bam_cleanup"),
    }

    if max_threads > 1:
        fixed_options["--max-threads"] = max_threads

    if paired_end:
        fixed_options["--paired-end"] = None

    if alt_aware or alt_optimize:
        convert.add_extra_files([InputFile(f"{in_reference}.alt")])

    if alt_optimize:
        fixed_options["--alt-optimize"] = None
        convert.add_extra_files([Executable("bwa-postalt.js")])

    convert.merge_options(
        user_options=options,
        fixed_options=fixed_options,
    )

    return convert


def _new_bwa_command(
    call: Iterable[AtomicFileTypes | str | int | Path],
    *,
    reference: PathTypes,
    index_ext: Iterable[str] = BWA_INDEX_EXT,
    index_iotype: type[InputFile | OutputFile] = InputFile,
    requirements: Iterable[versions.Requirement] = (BWA_VERSION,),
    stdout: int | str | Path | OutputFile | None = None,
) -> AtomicCmd:
    return AtomicCmd(
        call,
        extra_files=[index_iotype(os.fspath(reference) + ext) for ext in index_ext],
        requirements=requirements,
        stdout=stdout,
    )


@functools.lru_cache
def get_max_threads(reference: str, threads: int) -> int:
    """Returns the maximum number of threads to use when mapping against a
    given reference sequence. This is done since very little gain is obtained
    when using multiple threads for a small genome (e.g. < 1MB). If the
    reference falls below this size, only 1 thread is used (returned),
    otherwise the requested number of threads is returned.
    """
    if os.path.exists(reference) and os.path.getsize(reference) < 2**20:
        return 1

    return threads


def get_node_description(
    name: str,
    input_files_1: str,
    reference: str,
    input_files_2: str | None = None,
) -> str:
    reference = os.path.basename(reference)
    if reference.endswith((".fasta", ".fa")):
        reference = reference.rsplit(".", 1)[0]

    input_files = describe_paired_files(input_files_1, input_files_2 or ())

    return f"aligning {input_files} onto {reference} using {name}"
