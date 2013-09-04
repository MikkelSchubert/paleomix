#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
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

import pypeline.common.versions as versions

from pypeline.common.fileutils import \
    describe_files

from pypeline.node import \
     CommandNode
from pypeline.atomiccmd.command import \
     AtomicCmd
from pypeline.atomiccmd.sets import \
     ParallelCmds
from pypeline.nodes.picard import \
     concatenate_input_bams


# Number of reads to sample when running mapDamage
# FIXME: Remove and make node customizable instead
_MAPDAMAGE_MAX_READS = 100000

MAPDAMAGE_VERSION = versions.Requirement(call   = ("mapDamage", "--version"),
                                         search = r"(\d+)\.(\d+).(\d+)",
                                         checks = versions.GE(2, 0, 1))


class MapDamagePlotNode(CommandNode):
    def __init__(self, config, reference, input_files, output_directory, dependencies):
        cat_cmds, cat_obj = concatenate_input_bams(config, input_files)
        cmd_map = AtomicCmd(["mapDamage", "--no-stats",
                            "-n", _MAPDAMAGE_MAX_READS,
                             "-i", "-",
                             "-d", "%(TEMP_DIR)s",
                             "-r", "%(IN_REFERENCE)s"],
                             IN_STDIN        = cat_obj,
                             IN_REFERENCE    = reference,
                             OUT_FREQ_3p     = os.path.join(output_directory, "3pGtoA_freq.txt"),
                             OUT_FREQ_5p     = os.path.join(output_directory, "5pCtoT_freq.txt"),
                             OUT_COMP_USER   = os.path.join(output_directory, "dnacomp.txt"),
                             OUT_PLOT_FRAG   = os.path.join(output_directory, "Fragmisincorporation_plot.pdf"),
                             OUT_PLOT_LEN    = os.path.join(output_directory, "Length_plot.pdf"),
                             OUT_LENGTH      = os.path.join(output_directory, "lgdistribution.txt"),
                             OUT_MISINCORP   = os.path.join(output_directory, "misincorporation.txt"),
                             OUT_LOG         = os.path.join(output_directory, "Runtime_log.txt"),
                             CHECK_VERSION   = MAPDAMAGE_VERSION)

        description =  "<mapDamage (plots): %s -> '%s'>" \
          % (describe_files(input_files), output_directory)
        CommandNode.__init__(self,
                             command      = ParallelCmds(cat_cmds + [cmd_map]),
                             description  = description,
                             dependencies = dependencies)


class MapDamageModelNode(CommandNode):
    def __init__(self, reference, directory, dependencies):
        self._directory = directory

        command = AtomicCmd(["mapDamage", "--stats-only",
                             "-r", "%(IN_REFERENCE)s",
                             "-d", "%(TEMP_DIR)s"],
                             IN_REFERENCE         = reference,
                             TEMP_OUT_FREQ_3p     = "3pGtoA_freq.txt",
                             TEMP_OUT_FREQ_5p     = "5pCtoT_freq.txt",
                             TEMP_OUT_COMP_USER   = "dnacomp.txt",
                             TEMP_OUT_MISINCORP   = "misincorporation.txt",
                             TEMP_OUT_LOG         = "Runtime_log.txt",
                             OUT_COMP_GENOME      = os.path.join(directory, "dnacomp_genome.csv"),
                             OUT_MCMC_PROBS       = os.path.join(directory, "Stats_out_MCMC_correct_prob.csv"),
                             OUT_MCMC_HIST        = os.path.join(directory, "Stats_out_MCMC_hist.pdf"),
                             OUT_MCMC_ITER        = os.path.join(directory, "Stats_out_MCMC_iter.csv"),
                             OUT_MCMC_ITER_SUM    = os.path.join(directory, "Stats_out_MCMC_iter_summ_stat.csv"),
                             OUT_MCMC_POST_PRED   = os.path.join(directory, "Stats_out_MCMC_post_pred.pdf"),
                             OUT_MCMC_TRACE       = os.path.join(directory, "Stats_out_MCMC_trace.pdf"),
                             CHECK_VERSION        = MAPDAMAGE_VERSION)

        description =  "<mapDamage (model): %r>" % (directory,)
        CommandNode.__init__(self,
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)
        for fname in ("3pGtoA_freq.txt", "5pCtoT_freq.txt", "dnacomp.txt", "misincorporation.txt"):
            relpath = os.path.join(self._directory, fname)
            abspath = os.path.abspath(relpath)
            os.symlink(abspath, os.path.join(temp, fname))



class MapDamageRescaleNode(CommandNode):
    def __init__(self, config, reference, input_files, output_file, directory, dependencies):
        self._directory = directory

        cat_cmds, cat_obj = concatenate_input_bams(config, input_files)
        cmd_map = AtomicCmd(["mapDamage", "--rescale-only",
                             "-i", "-",
                             "-d", "%(TEMP_DIR)s",
                             "-r", "%(IN_REFERENCE)s",
                             "--rescale-out", "%(OUT_BAM)s"],
                             IN_STDIN        = cat_obj,
                             IN_REFERENCE    = reference,
                             TEMP_OUT_LOG    = "Runtime_log.txt",
                             TEMP_OUT_CSV    = "Stats_out_MCMC_correct_prob.csv",
                             OUT_BAM         = output_file,
                             CHECK_VERSION   = MAPDAMAGE_VERSION)

        description =  "<mapDamage (rescale): %s -> %r>" \
          % (describe_files(input_files), output_file)
        CommandNode.__init__(self,
                             command      = ParallelCmds(cat_cmds + [cmd_map]),
                             description  = description,
                             dependencies = dependencies)

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)
        for fname in ("Stats_out_MCMC_correct_prob.csv", ):
            relpath = os.path.join(self._directory, fname)
            abspath = os.path.abspath(relpath)
            os.symlink(abspath, os.path.join(temp, fname))
