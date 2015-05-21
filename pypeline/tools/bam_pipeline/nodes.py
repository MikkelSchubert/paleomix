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

from pypeline.node import \
    MetaNode
from pypeline.atomiccmd.command import \
    AtomicCmd
from pypeline.atomiccmd.sets import \
    ParallelCmds
from pypeline.atomiccmd.builder import \
    AtomicJavaCmdBuilder

from pypeline.nodes.picard import \
    PicardNode, \
    ValidateBAMNode, \
    MultiBAMInput, \
    MultiBAMInputNode
from pypeline.nodes.samtools import \
    BAMIndexNode
from pypeline.common.fileutils import \
    describe_files


def index_and_validate_bam(config, prefix, node, log_file=None):
    input_file, has_index = _get_input_file(node)
    if not has_index:
        node = BAMIndexNode(infile=input_file,
                            dependencies=node)

    validation_params = ValidateBAMNode.customize(config=config,
                                                  input_bam=input_file,
                                                  output_log=log_file,
                                                  dependencies=node)
    # Check MD tags against reference sequence
    # FIXME: Disabled due to issues with Picard/Samtools disagreeing,
    #   backwards compatibility. See the discussion at
    #     http://sourceforge.net/mailarchive/message.php?msg_id=31348639
    # validation_params.command.set_kwargs(IN_REF=prefix["Reference"])
    # validation_params.command.add_option("R", "%(IN_REF)s", sep="=")

    # Ignored since we may filter out misses and low-quality hits during
    # mapping, which leads to a large proportion of missing PE mates.
    validation_params.command.add_option("IGNORE", "MATE_NOT_FOUND",
                                         sep="=")
    # Ignored due to high rate of false positives for lanes with few hits,
    # where high-quality reads may cause mis-identification of qualities
    validation_params.command.add_option("IGNORE",
                                         "INVALID_QUALITY_FORMAT", sep="=")

    node = validation_params.build_node()
    return node


def _get_input_file(node):
    if isinstance(node, MetaNode):
        for subnode in node.subnodes:
            input_filename, has_index = _get_input_file(subnode)
            if input_filename:
                return input_filename, has_index

    input_filename, has_index = None, False
    for filename in node.output_files:
        if filename.lower().endswith(".bai"):
            has_index = True
        elif filename.lower().endswith(".bam"):
            input_filename = filename

    return input_filename, has_index


class CleanupBAMNode(PicardNode):
    def __init__(self, config, reference, input_bam, output_bam, tags,
                 min_mapq=0, filter_unmapped=False, dependencies=()):
        call = ["samtools", "view", "-bu"]
        if min_mapq > 0:
            call.append("-q%i" % min_mapq)
        if filter_unmapped:
            call.append("-F0x4")
        call.append("%(IN_BAM)s")

        flt = AtomicCmd(call,
                        IN_BAM=input_bam,
                        OUT_STDOUT=AtomicCmd.PIPE)

        jar_file = os.path.join(config.jar_root, "AddOrReplaceReadGroups.jar")
        params = AtomicJavaCmdBuilder(jar=jar_file,
                                      jre_options=config.jre_options)
        params.set_option("INPUT", "/dev/stdin", sep="=")
        params.set_option("OUTPUT", "%(TEMP_OUT_BAM)s", sep="=")
        params.set_option("COMPRESSION_LEVEL", "0", sep="=")
        params.set_option("SORT_ORDER", "coordinate", sep="=")

        for tag in ("ID", "SM", "LB", "PU", "PL"):
            params.set_option(tag, tags[tag], sep="=")

        params.set_kwargs(IN_STDIN=flt,
                          TEMP_OUT_BAM="bam.pipe")
        annotate = params.finalize()

        calmd = AtomicCmd(["samtools", "calmd", "-b",
                           "%(TEMP_IN_BAM)s", "%(IN_REF)s"],
                          IN_REF=reference,
                          TEMP_IN_BAM="bam.pipe",
                          OUT_STDOUT=output_bam)

        description = "<Cleanup BAM: %s -> '%s'>" \
            % (input_bam, output_bam)
        PicardNode.__init__(self,
                            command=ParallelCmds([flt, annotate, calmd]),
                            description=description,
                            dependencies=dependencies)

    def _setup(self, config, temp_root):
        PicardNode._setup(self, config, temp_root)
        os.mkfifo(os.path.join(temp_root, "bam.pipe"))
