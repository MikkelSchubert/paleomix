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

import paleomix.nodes.picard as picard

from paleomix.atomiccmd.command import \
    AtomicCmd
from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder
from paleomix.atomiccmd.sets import \
    ParallelCmds

from paleomix.nodes.picard import \
    PicardNode, \
    ValidateBAMNode
from paleomix.nodes.samtools import \
    BAMIndexNode


def index_and_validate_bam(config, prefix, node, log_file=None,
                           create_index=True):
    input_file, index_file = _get_input_files(node, prefix['IndexFormat'])
    if not index_file and create_index:
        node = BAMIndexNode(infile=input_file,
                            index_format=prefix['IndexFormat'],
                            dependencies=node)
        index_file, = node.output_files

    ignored_checks = [
        # Ignored since we may filter out misses and low-quality hits during
        # mapping, which leads to a large proportion of missing PE mates.
        "MATE_NOT_FOUND",
        # Ignored due to high rate of false positives for lanes with few hits,
        # where high-quality reads may cause mis-identification of qualities
        "INVALID_QUALITY_FORMAT"]

    if prefix['IndexFormat'] == '.csi':
        # CSI uses a different method for assigning BINs to records, which
        # Picard currently does not support.
        ignored_checks.append("INVALID_INDEXING_BIN")

    return ValidateBAMNode(config=config,
                           input_bam=input_file,
                           input_index=index_file,
                           ignored_checks=ignored_checks,
                           output_log=log_file,
                           dependencies=node)


def _get_input_files(node, index_format):
    index_filename = None
    input_filename = None
    for filename in node.output_files:
        if filename.lower().endswith(index_format):
            index_filename = True
        elif filename.lower().endswith(".bam"):
            input_filename = filename

    return input_filename, index_filename


class CleanupBAMNode(PicardNode):
    def __init__(self, config, reference, input_bam, output_bam, tags,
                 min_mapq=0, filter_unmapped=False, dependencies=()):
        flt_params = AtomicCmdBuilder(("samtools", "view", "-bu"),
                                      IN_BAM=input_bam,
                                      OUT_STDOUT=AtomicCmd.PIPE)

        if min_mapq:
            flt_params.set_option("-q", min_mapq, sep="")
        if filter_unmapped:
            flt_params.set_option("-F", "0x4", sep="")

        flt_params.add_value("%(IN_BAM)s")

        jar_params = picard.picard_command(config, "AddOrReplaceReadGroups")
        jar_params.set_option("INPUT", "/dev/stdin", sep="=")
        # Output is written to a named pipe, since the JVM may, in some cases,
        # emit warning messages to stdout, resulting in a malformed BAM.
        jar_params.set_option("OUTPUT", "%(TEMP_OUT_BAM)s", sep="=")
        jar_params.set_option("COMPRESSION_LEVEL", "0", sep="=")
        # Ensure that the BAM is sorted; this is required by the pipeline, and
        # needs to be done before calling calmd (avoiding pathologic runtimes).
        jar_params.set_option("SORT_ORDER", "coordinate", sep="=")

        # All tags are overwritten; ID is set since the default (e.g. '1')
        # causes problems with pysam due to type inference (is read as a length
        # 1 string, but written as a character).
        for tag in ("ID", "SM", "LB", "PU", "PL"):
            jar_params.set_option(tag, tags[tag], sep="=")

        jar_params.set_kwargs(IN_STDIN=flt_params,
                              TEMP_OUT_BAM="bam.pipe")

        calmd = AtomicCmdBuilder(["samtools", "calmd", "-b",
                                  "%(TEMP_IN_BAM)s", "%(IN_REF)s"],
                                 IN_REF=reference,
                                 TEMP_IN_BAM="bam.pipe",
                                 OUT_STDOUT=output_bam)

        commands = [cmd.finalize() for cmd in (flt_params, jar_params, calmd)]
        description = "<Cleanup BAM: %s -> '%s'>" \
            % (input_bam, output_bam)
        PicardNode.__init__(self,
                            command=ParallelCmds(commands),
                            description=description,
                            dependencies=dependencies)

    def _setup(self, config, temp_root):
        PicardNode._setup(self, config, temp_root)
        os.mkfifo(os.path.join(temp_root, "bam.pipe"))
