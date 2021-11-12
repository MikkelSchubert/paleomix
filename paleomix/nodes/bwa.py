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
import functools
import os

import paleomix.common.versions as versions
import paleomix.tools.factory as factory
from paleomix.common.command import (
    AtomicCmd,
    Executable,
    InputFile,
    OutputFile,
    ParallelCmds,
    TempOutputFile,
)
from paleomix.common.fileutils import describe_paired_files
from paleomix.node import CommandNode
from paleomix.nodes.samtools import SAMTOOLS_VERSION

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
    def __init__(self, input_file, dependencies=()):
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
            description="creating BWA index for %s" % (input_file,),
            dependencies=dependencies,
        )


class BWAMem2IndexNode(CommandNode):
    def __init__(self, input_file, dependencies=()):
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
            description="creating BWA MEM2 index for %s" % (input_file,),
            dependencies=dependencies,
        )


class BWABacktrack(CommandNode):
    def __init__(
        self,
        input_file,
        output_file,
        reference,
        threads=1,
        mapping_options={},
        dependencies=(),
    ):
        threads = _get_max_threads(reference, threads)
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
            description=_get_node_description(
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
        input_file_fq,
        input_file_sai,
        output_file,
        reference,
        threads=1,
        mapping_options={},
        cleanup_options={},
        dependencies=(),
    ):
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

        cleanup = _new_cleanup_command(
            stdin=samse,
            in_reference=reference,
            out_bam=output_file,
            max_threads=threads,
            options=cleanup_options,
        )

        CommandNode.__init__(
            self,
            command=ParallelCmds([samse, cleanup]),
            description=_get_node_description(
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
        input_file_fq_1,
        input_file_fq_2,
        input_file_sai_1,
        input_file_sai_2,
        output_file,
        reference,
        threads=1,
        mapping_options={},
        cleanup_options={},
        dependencies=(),
    ):
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

        cleanup = _new_cleanup_command(
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
            description=_get_node_description(
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
        input_file_1,
        output_file,
        reference,
        input_file_2=None,
        threads=1,
        algorithm="mem",
        alt_aware=False,
        alt_optimize=False,
        mapping_options={},
        cleanup_options={},
        dependencies=(),
    ):
        if algorithm in ("mem", "bwasw"):
            aln_call = ("bwa", algorithm, reference, InputFile(input_file_1))
            index_ext = BWA_INDEX_EXT
            requirements = [BWA_VERSION]
        elif algorithm == "mem2":
            aln_call = ("bwa-mem2", "mem", reference, InputFile(input_file_1))
            index_ext = BWA_MEM2_INDEX_EXT
            requirements = [BWA_MEM2_VERSION]
        else:
            raise NotImplementedError("BWA algorithm %r not implemented" % (algorithm,))

        aln = _new_bwa_command(
            aln_call,
            reference=reference,
            index_ext=index_ext,
            requirements=requirements,
            stdout=AtomicCmd.PIPE,
        )

        if input_file_2:
            aln.append(InputFile(input_file_2))

        if alt_aware or alt_optimize:
            if algorithm == "bwasw":
                raise NotImplementedError("bwasw is not ALT aware")

            aln.add_extra_files([InputFile(reference + ".alt")])

        threads = _get_max_threads(reference, threads)

        aln.merge_options(
            user_options=mapping_options,
            fixed_options={
                "-t": threads,
                # Mark alternative hits as secondary; required by e.g. Picard
                "-M": None,
            },
        )

        cleanup = _new_cleanup_command(
            stdin=aln,
            in_reference=reference,
            out_bam=output_file,
            max_threads=threads,
            paired_end=input_file_1 and input_file_2,
            alt_aware=alt_aware,
            alt_optimize=alt_optimize,
            options=cleanup_options,
        )

        self.out_bam = output_file

        CommandNode.__init__(
            self,
            command=ParallelCmds([aln, cleanup]),
            description=_get_node_description(
                name="BWA {}".format(algorithm),
                input_files_1=input_file_1,
                input_files_2=input_file_2,
                reference=reference,
            ),
            threads=threads,
            dependencies=dependencies,
        )


def _new_cleanup_command(
    stdin,
    in_reference,
    out_bam,
    max_threads=1,
    paired_end=False,
    alt_aware=False,
    alt_optimize=False,
    options={},
):
    convert = factory.new(
        "cleanup",
        stdin=stdin,
        stdout=out_bam,
        requirements=[
            SAMTOOLS_VERSION,
        ],
    )

    fixed_options = {
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
    call,
    reference,
    index_ext=BWA_INDEX_EXT,
    index_iotype=InputFile,
    requirements=(BWA_VERSION,),
    stdout=None,
):
    return AtomicCmd(
        call,
        extra_files=[index_iotype(reference + postfix) for postfix in index_ext],
        requirements=requirements,
        stdout=stdout,
    )


@functools.lru_cache()
def _get_max_threads(reference, threads):
    """Returns the maximum number of threads to use when mapping against a
    given reference sequence. This is done since very little gain is obtained
    when using multiple threads for a small genome (e.g. < 1MB). If the
    reference falls below this size, only 1 thread is used (returned),
    otherwise the requested number of threads is returned.
    """
    if os.path.exists(reference) and os.path.getsize(reference) < 2 ** 20:
        return 1

    return threads


def _get_node_description(name, input_files_1, reference, input_files_2=None):
    reference = os.path.basename(reference)
    if reference.endswith(".fasta") or reference.endswith(".fa"):
        reference = reference.rsplit(".", 1)[0]

    input_files = describe_paired_files(input_files_1, input_files_2 or ())

    return "aligning {} onto {} using {}".format(input_files, reference, name)
