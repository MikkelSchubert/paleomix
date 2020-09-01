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

from paleomix.atomiccmd.builder import (
    AtomicCmdBuilder,
    apply_options,
)
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import ParallelCmds
from paleomix.common.fileutils import describe_paired_files
from paleomix.node import CommandNode, NodeError
from paleomix.nodes.samtools import SAMTOOLS_VERSION


BWA_VERSION = versions.Requirement(
    name="BWA",
    call=("bwa",),
    search=r"Version: (\d+)\.(\d+)\.(\d+)",
    checks=versions.GE(0, 7, 15),
)


class BWAIndexNode(CommandNode):
    def __init__(self, input_file, prefix=None, dependencies=()):
        prefix = prefix if prefix else input_file
        builder = _new_bwa_command(
            ("bwa", "index", "%(IN_FILE)s", "-p", "%(TEMP_OUT_PREFIX)s"),
            prefix,
            iotype="OUT",
            IN_FILE=input_file,
            TEMP_OUT_PREFIX=os.path.basename(prefix),
        )

        CommandNode.__init__(
            self,
            command=builder.finalize(),
            description="creating BWA index for %s" % (input_file,),
            dependencies=dependencies,
        )


class BWABacktrack(CommandNode):
    def __init__(
        self,
        input_file,
        output_file,
        reference,
        prefix,
        threads=1,
        mapping_options={},
        dependencies=(),
    ):
        threads = _get_max_threads(reference, threads)

        aln = _new_bwa_command(
            ("bwa", "aln"), prefix, IN_FILE=input_file, OUT_STDOUT=output_file,
        )
        aln.add_value(prefix)
        aln.add_value("%(IN_FILE)s")
        aln.set_option("-t", threads)

        apply_options(aln, mapping_options)

        description = _get_node_description(
            name="BWA",
            algorithm="Backtrack",
            input_files_1=input_file,
            prefix=prefix,
            threads=threads,
        )

        CommandNode.__init__(
            self,
            command=aln.finalize(),
            description=description,
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
        prefix,
        mapping_options={},
        cleanup_options={},
        dependencies=(),
    ):
        samse = _new_bwa_command(
            ("bwa", "samse"),
            prefix,
            IN_FILE_SAI=input_file_sai,
            IN_FILE_FQ=input_file_fq,
            OUT_STDOUT=AtomicCmd.PIPE,
        )
        samse.add_value(prefix)
        samse.add_value("%(IN_FILE_SAI)s")
        samse.add_value("%(IN_FILE_FQ)s")

        cleanup = _new_cleanup_command(samse, output_file, reference)

        apply_options(samse, mapping_options)
        apply_options(cleanup, cleanup_options)

        CommandNode.__init__(
            self,
            command=ParallelCmds([samse.finalize(), cleanup.finalize()]),
            description=_get_node_description(
                name="BWA",
                algorithm="samse",
                input_files_1=input_file_fq,
                prefix=prefix,
            ),
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
        prefix,
        mapping_options={},
        cleanup_options={},
        dependencies=(),
    ):
        sampe = _new_bwa_command(
            (
                "bwa",
                "sampe",
                prefix,
                "%(IN_SAI_1)s",
                "%(IN_SAI_2)s",
                "%(IN_FQ_1)s",
                "%(IN_FQ_2)s",
            ),
            prefix,
            IN_SAI_1=input_file_sai_1,
            IN_SAI_2=input_file_sai_2,
            IN_FQ_1=input_file_fq_1,
            IN_FQ_2=input_file_fq_2,
            OUT_STDOUT=AtomicCmd.PIPE,
        )

        cleanup = _new_cleanup_command(sampe, output_file, reference, paired_end=True)

        apply_options(sampe, mapping_options)
        apply_options(cleanup, cleanup_options)

        CommandNode.__init__(
            self,
            command=ParallelCmds([sampe.finalize(), cleanup.finalize()]),
            description=_get_node_description(
                name="BWA",
                algorithm="sampe",
                input_files_1=input_file_fq_1,
                input_files_2=input_file_fq_2,
                prefix=prefix,
            ),
            dependencies=dependencies,
        )


class BWAAlgorithmNode(CommandNode):
    def __init__(
        self,
        input_file_1,
        output_file,
        reference,
        prefix,
        input_file_2=None,
        threads=1,
        algorithm="mem",
        mapping_options={},
        cleanup_options={},
        dependencies=(),
    ):
        if algorithm not in ("mem", "bwasw"):
            raise NotImplementedError("BWA algorithm %r not implemented" % (algorithm,))

        threads = _get_max_threads(reference, threads)

        aln = _new_bwa_command(
            ("bwa", algorithm, prefix, "%(IN_FILE_1)s"),
            prefix,
            IN_FILE_1=input_file_1,
            OUT_STDOUT=AtomicCmd.PIPE,
        )

        if input_file_2:
            aln.add_value("%(IN_FILE_2)s")
            aln.set_kwargs(IN_FILE_2=input_file_2)

        aln.set_option("-t", threads)
        # Mark alternative hits as secondary; required by e.g. Picard
        aln.set_option("-M")

        cleanup = _new_cleanup_command(
            aln, output_file, reference, paired_end=input_file_1 and input_file_2
        )

        apply_options(aln, mapping_options)
        apply_options(cleanup, cleanup_options)

        description = _get_node_description(
            name="BWA",
            algorithm="%s (%s)" % (algorithm, "PE" if input_file_2 else "SE"),
            input_files_1=input_file_1,
            input_files_2=input_file_2,
            prefix=prefix,
            threads=threads,
        )

        CommandNode.__init__(
            self,
            command=ParallelCmds([aln.finalize(), cleanup.finalize()]),
            description=description,
            threads=threads,
            dependencies=dependencies,
        )


def _new_cleanup_command(stdin, output_file, reference, paired_end=False):
    convert = factory.new("cleanup")
    convert.set_option("--fasta", "%(IN_FASTA_REF)s")
    convert.set_option("--temp-prefix", "%(TEMP_OUT_PREFIX)s")
    convert.set_kwargs(
        IN_STDIN=stdin,
        IN_FASTA_REF=reference,
        OUT_STDOUT=output_file,
        TEMP_OUT_PREFIX="bam_cleanup",
        CHECK_SAMTOOLS=SAMTOOLS_VERSION,
    )

    if paired_end:
        convert.set_option("--paired-end")

    return convert


def _new_bwa_command(call, prefix, iotype="IN", **kwargs):
    _check_bwa_prefix(prefix)

    kwargs["CHECK_BWA"] = BWA_VERSION
    for postfix in ("amb", "ann", "bwt", "pac", "sa"):
        kwargs["%s_PREFIX_%s" % (iotype, postfix.upper())] = prefix + "." + postfix

    return AtomicCmdBuilder(call, **kwargs)


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


@functools.lru_cache()
def _check_bwa_prefix(prefix):
    """Checks that a given prefix is compatible with the currently required version of
    BWA. Older index files are incompatible with BWA v0.7.x, but may be identified by
    the presense of a small number of additional files not present when an index is
    produced using BWA v0.7.x.
    """
    if any(os.path.exists(prefix + ext) for ext in (".rbwt", ".rpac", ".rsa.")):
        filenames = "\n".join(
            "    %s.%s" % (prefix, ext)
            for ext in ("amb", "ann", "bwt", "pac", "sa", "rbwt", "rpac", "rsa")
        )

        raise NodeError(
            "Prefix appears to be created using BWA v0.5.x or older, but PALEOMIX only "
            "supports BWA v0.7.x or later.\nPlease remove the following files to allow "
            "PALEOMIX to re-index the FASTA file:\n%s" % (filenames,)
        )


def _get_node_description(
    name, input_files_1, input_files_2=None, algorithm=None, prefix=None, threads=1
):
    prefix = os.path.basename(prefix)
    if prefix.endswith(".fasta") or prefix.endswith(".fa"):
        prefix = prefix.rsplit(".", 1)[0]

    input_files = describe_paired_files(input_files_1, input_files_2 or ())
    info = ["alignment of", input_files, "onto", prefix, "using", name]
    if algorithm is not None:
        info.append(algorithm)

    if threads > 1:
        info.append("using %i threads" % (threads,))

    return " ".join(info)
