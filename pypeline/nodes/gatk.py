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

import pypeline.common.fileutils as \
    fileutils
from pypeline.node import \
    MetaNode, \
    CommandNode
from pypeline.atomiccmd.command import \
    AtomicCmd
from pypeline.atomiccmd.builder import \
    AtomicJavaCmdBuilder
from pypeline.atomiccmd.sets import \
    ParallelCmds
from pypeline.common.fileutils import \
    swap_ext, \
    describe_files
from pypeline.common.utilities import \
    safe_coerce_to_tuple
from pypeline.nodes.samtools import \
    SAMTOOLS_VERSION

import pypeline.common.versions as versions


def _get_gatk_version_check(config):
    """Returns a version-check object for the "GenomeAnalysisTK.jar" located at
    config.jar_root; for now, this check only serves to verify that the JAR can
    be executed, which may not be the case if the JRE is outdated.
    """
    jar_file = os.path.join(config.jar_root, "GenomeAnalysisTK.jar")
    if jar_file not in _GATK_VERSION:
        params = AtomicJavaCmdBuilder(jar_file,
                                      temp_root=config.temp_root,
                                      jre_options=config.jre_options)
        params.add_value("--version")

        # Any version is fine; for now just catch old JREs
        requirement = versions.Requirement(call=params.finalized_call,
                                           name="GenomeAnalysisTK",
                                           search=r"^(\d+)\.(\d+)",
                                           checks=versions.Any())
        _GATK_VERSION[jar_file] = requirement
    return _GATK_VERSION[jar_file]
_GATK_VERSION = {}


class _IndelTrainerNode(CommandNode):
    def __init__(self, config, reference, infiles, outfile, dependencies=()):
        infiles = safe_coerce_to_tuple(infiles)
        jar_file = os.path.join(config.jar_root, "GenomeAnalysisTK.jar")
        command = AtomicJavaCmdBuilder(jar_file,
                                       jre_options=config.jre_options)
        command.set_option("-T", "RealignerTargetCreator")
        command.set_option("-R", "%(IN_REFERENCE)s")
        command.set_option("-o", "%(OUT_INTERVALS)s")

        _set_input_files(command, infiles)
        command.set_kwargs(IN_REFERENCE=reference,
                           IN_REF_DICT=fileutils.swap_ext(reference, ".dict"),
                           OUT_INTERVALS=outfile,
                           CHECK_GATK=_get_gatk_version_check(config))

        description = "<Indel Realigner (training): %s -> %r>" \
            % (describe_files(infiles), outfile)
        CommandNode.__init__(self,
                             description=description,
                             command=command.finalize(),
                             dependencies=dependencies)


class _IndelRealignerNode(CommandNode):
    def __init__(self, config, reference, intervals, infiles, outfile,
                 dependencies=()):
        self._basename = os.path.basename(outfile)

        infiles = safe_coerce_to_tuple(infiles)
        jar_file = os.path.join(config.jar_root, "GenomeAnalysisTK.jar")
        command = AtomicJavaCmdBuilder(jar_file,
                                       jre_options=config.jre_options)
        command.set_option("-T", "IndelRealigner")
        command.set_option("-R", "%(IN_REFERENCE)s")
        command.set_option("-targetIntervals", "%(IN_INTERVALS)s")
        command.set_option("-o", "%(OUT_BAMFILE)s")
        command.set_option("--bam_compression", 0)
        command.set_option("--disable_bam_indexing")
        _set_input_files(command, infiles)

        command.set_kwargs(IN_REFERENCE=reference,
                           IN_REF_DICT=fileutils.swap_ext(reference, ".dict"),
                           IN_INTERVALS=intervals,
                           OUT_BAMFILE=outfile,
                           CHECK_GATK=_get_gatk_version_check(config))

        calmd = AtomicCmd(["samtools", "calmd", "-b",
                           "%(TEMP_IN_BAM)s", "%(IN_REF)s"],
                          TEMP_IN_BAM=self._basename,
                          IN_REF=reference,
                          TEMP_OUT_STDOUT=self._basename + ".calmd",
                          CHECK_VERSION=SAMTOOLS_VERSION)

        description = "<Indel Realigner (aligning): %s -> %r>" \
            % (describe_files(infiles), outfile)
        CommandNode.__init__(self,
                             description=description,
                             command=ParallelCmds([command.finalize(), calmd]),
                             dependencies=dependencies)

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)
        os.mkfifo(os.path.join(temp, self._basename))

    def _teardown(self, config, temp):
        os.rename(os.path.join(temp, self._basename) + ".calmd",
                  os.path.join(temp, self._basename))

        CommandNode._teardown(self, config, temp)


class IndelRealignerNode(MetaNode):
    def __init__(self, config, reference, infiles, outfile, intervals=None,
                 dependencies=()):
        if not intervals:
            intervals = outfile + ".intervals"

        infiles = safe_coerce_to_tuple(infiles)
        trainer = _IndelTrainerNode(config=config,
                                    reference=reference,
                                    infiles=infiles,
                                    outfile=intervals,
                                    dependencies=dependencies)
        aligner = _IndelRealignerNode(config=config,
                                      reference=reference,
                                      intervals=intervals,
                                      infiles=infiles,
                                      outfile=outfile,
                                      dependencies=trainer)

        description = "<GATK Indel Realigner: %i files -> '%s'>" \
            % (len(infiles), outfile)

        MetaNode.__init__(self,
                          description=description,
                          subnodes=[trainer, aligner],
                          dependencies=dependencies)


def _set_input_files(command, input_files):
    keys = {}
    for (index, filename) in enumerate(input_files):
        command.add_option("-I", "%%(IN_BAMFILE_%02i)s" % index)
        keys["IN_BAMFILE_%02i" % index] = filename
        keys["IN_BAIFILE_%02i" % index] = swap_ext(filename, ".bai")

    command.set_kwargs(**keys)
