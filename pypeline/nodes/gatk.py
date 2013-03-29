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

import pypeline.common.fileutils as fileutils
from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds
from pypeline.common.fileutils import swap_ext
from pypeline.common.utilities import safe_coerce_to_tuple


class _IndelTrainerNode(CommandNode):
    def __init__(self, config, reference, infiles, outfile, dependencies = ()):
        jar   = os.path.join(config.jar_root, "GenomeAnalysisTK.jar")
        call  = ["java", "-jar", jar,
                 "-T", "RealignerTargetCreator",
                 "-R", "%(IN_REFERENCE)s",
                 "-o", "%(OUT_INTERVALS)s"]
        keys = _update_file_keys_and_call(call, infiles)

        command = AtomicCmd(call,
                            IN_REFERENCE  = reference,
                            IN_REF_DICT   = fileutils.swap_ext(reference, ".dict"),
                            OUT_INTERVALS = outfile,
                            AUX_JAR       = jar,
                            **keys)

        description = "<Train Indel Realigner: %i file(s) -> '%s'>" \
            % (len(infiles), outfile)

        CommandNode.__init__(self,
                             description = description,
                             command = command,
                             dependencies = dependencies)



class _IndelRealignerNode(CommandNode):
    def __init__(self, config, reference, intervals, infiles, outfile, dependencies = ()):
        self._basename = os.path.basename(outfile)

        jar   = os.path.join(config.jar_root, "GenomeAnalysisTK.jar")
        call  = ["java", "-jar", jar,
                 "-T", "IndelRealigner",
                 "-R", "%(IN_REFERENCE)s",
                 "-targetIntervals", "%(IN_INTERVALS)s",
                 "-o", "%(OUT_BAMFILE)s",
                 "--bam_compression", 0,
                 "--disable_bam_indexing"]
        keys = _update_file_keys_and_call(call, infiles)

        command = AtomicCmd(call,
                            IN_REFERENCE = reference,
                            IN_REF_DICT  = fileutils.swap_ext(reference, ".dict"),
                            IN_INTERVALS = intervals,
                            OUT_BAMFILE  = outfile,
                            AUX_JAR      = jar,
                            **keys)

        calmd   = AtomicCmd(["samtools", "calmd", "-b", "%(TEMP_IN_BAM)s", "%(IN_REF)s"],
                            TEMP_IN_BAM     = self._basename,
                            IN_REF          = reference,
                            TEMP_OUT_STDOUT = self._basename + ".calmd")

        description = "<Indel Realign: %i file(s) -> '%s'>" \
            % (len(infiles), outfile)

        CommandNode.__init__(self,
                             description  = description,
                             command      = ParallelCmds([command, calmd]),
                             dependencies = dependencies)


    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)
        os.mkfifo(os.path.join(temp, self._basename))


    def _teardown(self, config, temp):
        os.rename(os.path.join(temp, self._basename) + ".calmd",
                  os.path.join(temp, self._basename))

        CommandNode._teardown(self, config, temp)


class IndelRealignerNode(MetaNode):
    def __init__(self, config, reference, infiles, outfile, intervals = None, dependencies = ()):
        if not intervals:
            intervals = outfile + ".intervals"

        infiles = safe_coerce_to_tuple(infiles)
        trainer = _IndelTrainerNode(config         = config,
                                    reference      = reference,
                                    infiles        = infiles,
                                    outfile        = intervals,
                                    dependencies   = dependencies)
        aligner = _IndelRealignerNode(config       = config,
                                      reference    = reference,
                                      intervals    = intervals,
                                      infiles      = infiles,
                                      outfile      = outfile,
                                      dependencies = trainer)

        MetaNode.__init__(self,
                          description  = "<GATK Indel Realigner: %i files -> '%s'>"
                                             % (len(infiles), outfile),
                          subnodes     = [trainer, aligner],
                          dependencies = dependencies)


def _update_file_keys_and_call(call, input_files):
    keys = {}
    for (index, filename) in enumerate(input_files):
        call.extend(("-I", "%%(IN_BAMFILE_%02i)s" % index))
        keys["IN_BAMFILE_%02i" % index] = filename
        keys["IN_BAIFILE_%02i" % index] = swap_ext(filename, ".bai")

    return keys
