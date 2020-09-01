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

import paleomix.common.rtools as rtools
import paleomix.common.versions as versions

from paleomix.common.fileutils import describe_files

from paleomix.node import NodeError, CommandNode
from paleomix.nodes.samtools import merge_bam_files_command
from paleomix.atomiccmd.builder import AtomicCmdBuilder, apply_options
from paleomix.atomiccmd.sets import ParallelCmds


MAPDAMAGE_VERSION = versions.Requirement(
    call=("mapDamage", "--version"),
    search=r"(\d+)\.(\d+).(\d+)",
    checks=versions.GE(2, 2, 1),
)

RSCRIPT_VERSION = versions.Requirement(
    call=("Rscript", "--version"),
    search=r"version (\d+)\.(\d+).(\d+)",
    checks=versions.GE(3, 3, 3),
    priority=10,
)


class MapDamagePlotNode(CommandNode):
    def __init__(
        self,
        reference,
        input_files,
        output_directory,
        title="mapDamage",
        options={},
        dependencies=(),
    ):
        merge = merge_bam_files_command(input_files)
        command = AtomicCmdBuilder(
            [
                "mapDamage",
                "--no-stats",
                # Prevent references with many contigs from using excessive
                # amounts of memory, at the cost of per-contig statistics:
                "--merge-reference-sequences",
                "-t",
                title,
                "-i",
                "-",
                "-d",
                "%(TEMP_DIR)s",
                "-r",
                "%(IN_REFERENCE)s",
            ],
            IN_STDIN=merge,
            IN_REFERENCE=reference,
            OUT_FREQ_3p=os.path.join(output_directory, "3pGtoA_freq.txt"),
            OUT_FREQ_5p=os.path.join(output_directory, "5pCtoT_freq.txt"),
            OUT_COMP_USER=os.path.join(output_directory, "dnacomp.txt"),
            OUT_PLOT_FRAG=os.path.join(
                output_directory, "Fragmisincorporation_plot.pdf"
            ),
            OUT_PLOT_LEN=os.path.join(output_directory, "Length_plot.pdf"),
            OUT_LENGTH=os.path.join(output_directory, "lgdistribution.txt"),
            OUT_MISINCORP=os.path.join(output_directory, "misincorporation.txt"),
            OUT_LOG=os.path.join(output_directory, "Runtime_log.txt"),
            TEMP_OUT_STDOUT="pipe_mapDamage.stdout",
            TEMP_OUT_STDERR="pipe_mapDamage.stderr",
            CHECK_RSCRIPT=RSCRIPT_VERSION,
            CHECK_MAPDAMAGE=MAPDAMAGE_VERSION,
        )

        apply_options(command, options)

        CommandNode.__init__(
            self,
            command=ParallelCmds([merge, command.finalize()]),
            description="creating mapDamage plots from %s"
            % (describe_files(merge.input_files),),
            dependencies=dependencies,
        )

    def _teardown(self, config, temp):
        # No Length_plot.pdf file is written if there are no SE reads in the
        # input_file. In that case, we write a dummy PDF to ensure that all
        # expected files exist.
        err_message = "No length distributions are available"
        with open(os.path.join(temp, "pipe_mapDamage.stderr")) as in_handle:
            if any(line.startswith(err_message) for line in in_handle):
                fpath = os.path.join(temp, "Length_plot.pdf")
                with open(fpath, "w") as out_handle:
                    out_handle.write(_DUMMY_LENGTH_PLOT_PDF)

        CommandNode._teardown(self, config, temp)


class MapDamageModelNode(CommandNode):
    def __init__(self, reference, directory, options={}, dependencies=()):
        command = AtomicCmdBuilder(
            [
                "mapDamage",
                "--stats-only",
                "-r",
                "%(IN_REFERENCE)s",
                "-d",
                "%(TEMP_DIR)s",
            ],
            IN_REFERENCE=reference,
            TEMP_OUT_FREQ_3p="3pGtoA_freq.txt",
            TEMP_OUT_FREQ_5p="5pCtoT_freq.txt",
            TEMP_OUT_COMP_USER="dnacomp.txt",
            TEMP_OUT_MISINCORP="misincorporation.txt",
            TEMP_OUT_LOG="Runtime_log.txt",
            TEMP_OUT_STDOUT="pipe_mapDamage.stdout",
            TEMP_OUT_STDERR="pipe_mapDamage.stderr",
            OUT_COMP_GENOME=os.path.join(directory, "dnacomp_genome.csv"),
            OUT_MCMC_PROBS=os.path.join(directory, "Stats_out_MCMC_correct_prob.csv"),
            OUT_MCMC_HIST=os.path.join(directory, "Stats_out_MCMC_hist.pdf"),
            OUT_MCMC_ITER=os.path.join(directory, "Stats_out_MCMC_iter.csv"),
            OUT_MCMC_ITERSUM=os.path.join(
                directory, "Stats_out_MCMC_iter_summ_stat.csv"
            ),
            OUT_MCMC_POSTPRED=os.path.join(directory, "Stats_out_MCMC_post_pred.pdf"),
            OUT_MCMC_TRACE=os.path.join(directory, "Stats_out_MCMC_trace.pdf"),
            CHECK_RSCRIPT=RSCRIPT_VERSION,
            CHECK_MAPDAMAGE=MAPDAMAGE_VERSION,
            CHECK_R_INLINE=rtools.requirement("inline"),
            CHECK_R_GGPLOT2=rtools.requirement("ggplot2"),
            CHECK_R_RCPP=rtools.requirement("Rcpp"),
            CHECK_R_GAM=rtools.requirement("gam"),
            CHECK_R_RCPPGSL=rtools.requirement("RcppGSL"),
        )

        apply_options(command, options)

        self._directory = directory

        CommandNode.__init__(
            self,
            command=command.finalize(),
            description="creating mapDamage model from %r" % (directory,),
            dependencies=dependencies,
        )

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)
        for fname in (
            "3pGtoA_freq.txt",
            "5pCtoT_freq.txt",
            "dnacomp.txt",
            "misincorporation.txt",
        ):
            relpath = os.path.join(self._directory, fname)
            abspath = os.path.abspath(relpath)
            os.symlink(abspath, os.path.join(temp, fname))

    def _run(self, config, temp):
        try:
            CommandNode._run(self, config, temp)
        except NodeError as error:
            err_message = "DNA damage levels are too low"
            if self._command.join() == [1]:
                fpath = os.path.join(temp, "pipe_mapDamage.stdout")
                with open(fpath) as handle:
                    for line in handle:
                        if err_message in line:
                            line = line.strip().replace("Warning:", "ERROR:")
                            error = NodeError("%s\n\n%s" % (error, line))
                            break
            raise error


class MapDamageRescaleNode(CommandNode):
    def __init__(
        self,
        reference,
        input_files,
        output_file,
        directory,
        options={},
        dependencies=(),
    ):
        stats_out_fname = "Stats_out_MCMC_correct_prob.csv"

        merge = merge_bam_files_command(input_files)
        command = AtomicCmdBuilder(
            [
                "mapDamage",
                "--rescale-only",
                "-i",
                "-",
                "-d",
                "%(TEMP_DIR)s",
                "-r",
                "%(IN_REFERENCE)s",
                "--rescale-out",
                "%(OUT_BAM)s",
            ],
            IN_STDIN=merge,
            IN_REFERENCE=reference,
            TEMP_OUT_LOG="Runtime_log.txt",
            TEMP_OUT_CSV=stats_out_fname,
            OUT_BAM=output_file,
            CHECK_VERSION=MAPDAMAGE_VERSION,
        )

        apply_options(command, options)

        self._directory = directory

        CommandNode.__init__(
            self,
            command=ParallelCmds([merge, command.finalize()]),
            description="rescaling %s using mapDamage"
            % (describe_files(merge.input_files),),
            dependencies=dependencies,
        )

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)
        for fname in ("Stats_out_MCMC_correct_prob.csv",):
            relpath = os.path.join(self._directory, fname)
            abspath = os.path.abspath(relpath)
            os.symlink(abspath, os.path.join(temp, fname))


# Minimal PDF written if Length_plot.pdf wasn't generated
_DUMMY_LENGTH_PLOT_PDF = """%PDF-1.4

1 0 obj
 <</Type /Font /Subtype /Type1 /Encoding /WinAnsiEncoding /BaseFont /Courier >>
endobj

2 0 obj
 <</Parent 4 0 R /MediaBox[0 0 450 50] /Type /Page /Contents[3 0 R ] /Resources 5 0 R >>
endobj

3 0 obj
 <</Length 138 >>
stream
  BT
    /F0 18 Tf
    20 10 Td
    (Input file(s) did not contain SE reads.) Tj
    0 20 Td
    (Length_plot.pdf not generated:) Tj
  ET
endstream
endobj

4 0 obj
 <</Type /Pages /Count 1 /Kids[2 0 R ]>>
endobj

5 0 obj
 <</ProcSet[/PDF /Text] /Font <</F0 1 0 R >>
>>
endobj

6 0 obj
 <</Type /Catalog /Pages 4 0 R >>
endobj

xref
0 7
0000000000 65535 f
0000000010 00000 n
0000000106 00000 n
0000000211 00000 n
0000000400 00000 n
0000000457 00000 n
0000000521 00000 n
trailer
 <</Size 7 /Root 6 0 R >>

startxref
571
%%EOF
"""
