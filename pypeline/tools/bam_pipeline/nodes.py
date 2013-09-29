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

from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
from pypeline.atomiccmd.builder import AtomicJavaCmdBuilder

from pypeline.nodes.picard import ValidateBAMNode, concatenate_input_bams
from pypeline.nodes.samtools import BAMIndexNode
from pypeline.common.fileutils import describe_files




class FilterCollapsedBAMNode(CommandNode):
    def __init__(self, config, input_bams, output_bam, keep_dupes, dependencies = ()):
        cat_cmds, cat_obj = concatenate_input_bams(config, input_bams)
        rmdup_call = ["bam_rmdup_collapsed"]
        if not keep_dupes:
            rmdup_call.append("--remove-duplicates")

        filteruniq = AtomicCmd(rmdup_call,
                               IN_STDIN   = cat_obj,
                               OUT_STDOUT = output_bam)

        command     = ParallelCmds(cat_cmds + [filteruniq])
        description =  "<FilterCollapsedBAM: %s>" % (describe_files(input_bams),)
        CommandNode.__init__(self,
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


class IndexAndValidateBAMNode(MetaNode):
    def __init__(self, config, prefix, node, log_file = None):
        input_file, has_index = self._get_input_file(node)
        subnodes, dependencies = [node], node.dependencies
        if not has_index:
            node = BAMIndexNode(infile       = input_file,
                                dependencies = node)
            subnodes.append(node)

        validation_params = ValidateBAMNode.customize(config       = config,
                                                      input_bam    = input_file,
                                                      output_log   = log_file,
                                                      dependencies = node)
        # Check MD tags against reference sequence
        # FIXME: Disabled due to issues with Picard/Samtools disagreeing, backwards compatibility.
        #        See: http://sourceforge.net/mailarchive/message.php?msg_id=31348639
        #        validation_params.command.set_kwargs(IN_REFERENCE = prefix["Reference"])
        #        validation_params.command.add_option("R", "%(IN_REFERENCE)s", sep = "=")
        # Ignored since we filter out misses and low-quality hits during mapping, which
        # leads to a large proportion of missing mates for PE reads.
        validation_params.command.add_option("IGNORE", "MATE_NOT_FOUND", sep = "=")
        # Ignored due to high rate of false positives for lanes with few hits, where
        # high-quality reads may case ValidateSamFile to mis-identify the qualities
        validation_params.command.add_option("IGNORE", "INVALID_QUALITY_FORMAT", sep = "=")
        subnodes.append(validation_params.build_node())

        description = "<w/Validation: " + str(subnodes[0])[1:]
        MetaNode.__init__(self,
                          description  = description,
                          subnodes     = subnodes,
                          dependencies = dependencies)


    @classmethod
    def _get_input_file(cls, node):
        if isinstance(node, MetaNode):
            for subnode in node.subnodes:
                input_filename, has_index = cls._get_input_file(subnode)
                if input_filename:
                    return input_filename, has_index

        input_filename, has_index = None, False
        for filename in node.output_files:
            if filename.lower().endswith(".bai"):
                has_index = True
            elif filename.lower().endswith(".bam"):
                input_filename = filename

        return input_filename, has_index


class CleanupBAMNode(CommandNode):
    def __init__(self, config, reference, input_bam, output_bam, tags, min_mapq = 0, dependencies = ()):
        flt = AtomicCmd(["samtools", "view", "-bu", "-F0x4", "-q%i" % min_mapq, "%(IN_BAM)s"],
                        IN_BAM  = input_bam,
                        OUT_STDOUT = AtomicCmd.PIPE)

        jar_file = os.path.join(config.jar_root, "AddOrReplaceReadGroups.jar")
        params = AtomicJavaCmdBuilder(config, jar_file)
        params.set_option("INPUT", "/dev/stdin", sep = "=")
        params.set_option("OUTPUT", "/dev/stdout", sep = "=")
        params.set_option("QUIET", "true", sep = "=")
        params.set_option("COMPRESSION_LEVEL", "0", sep = "=")
        params.set_option("SORT_ORDER", "coordinate", sep = "=")

        for tag in ("SM", "LB", "PU", "PL"):
            params.set_option(tag, tags[tag], sep = "=")

        params.set_kwargs(IN_STDIN   = flt,
                         OUT_STDOUT = AtomicCmd.PIPE)
        annotate = params.finalize()

        calmd = AtomicCmd(["samtools", "calmd", "-b", "-", "%(IN_REF)s"],
                          IN_REF   = reference,
                          IN_STDIN = annotate,
                          OUT_STDOUT = output_bam)

        description =  "<Cleanup BAM: %s -> '%s'>" \
            % (input_bam, output_bam)
        CommandNode.__init__(self,
                             command      = ParallelCmds([flt, annotate, calmd]),
                             description  = description,
                             dependencies = dependencies)


