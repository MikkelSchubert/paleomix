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
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds
from pypeline.atomicparams import AtomicJavaParams

from pypeline.nodes.picard import ValidateBAMNode
from pypeline.nodes.samtools import BAMIndexNode
from pypeline.common.utilities import safe_coerce_to_tuple
from pypeline.common.fileutils import swap_ext, add_postfix


# Number of reads to sample when running mapDamage
_MAPDAMAGE_MAX_READS = 100000




class MapDamageNode(CommandNode):
    def __init__(self, reference, input_files, output_directory, dependencies):
        input_files = safe_coerce_to_tuple(input_files)

        # Workaround to simplify handling of 1/2+ files (samtools cat requires 2+!)
        cat_cmd   = ["cat"] if (len(input_files) == 1) else ["samtools", "cat"]
        cat_files = {"OUT_STDOUT" : AtomicCmd.PIPE}
        for (index, filename) in enumerate(safe_coerce_to_tuple(input_files)):
            cat_cmd.append("%%(IN_FILE_%02i)s" % (index,))
            cat_files["IN_FILE_%02i" % (index,)] = filename
        cmd_cat = AtomicCmd(cat_cmd, **cat_files)

        cmd_map = AtomicCmd(["mapDamage", "--no-stats",
                            "-n", _MAPDAMAGE_MAX_READS,
                             "-i", "-",
                             "-d", "%(TEMP_DIR)s",
                             "-r", reference],
                            IN_STDIN        = cmd_cat,
                            OUT_FREQ_3p     = os.path.join(output_directory, "3pGtoA_freq.txt"),
                            OUT_FREQ_5p     = os.path.join(output_directory, "5pCtoT_freq.txt"),
                            OUT_COMP_GENOME = os.path.join(output_directory, "dnacomp_genome.csv"),
                            OUT_COMP_USER   = os.path.join(output_directory, "dnacomp.txt"),
                            OUT_PLOT_FRAG   = os.path.join(output_directory, "Fragmisincorporation_plot.pdf"),
                            OUT_PLOT_LEN    = os.path.join(output_directory, "Length_plot.pdf"),
                            OUT_LENGTH      = os.path.join(output_directory, "lgdistribution.txt"),
                            OUT_MISINCORP   = os.path.join(output_directory, "misincorporation.txt"),
                            OUT_LOG         = os.path.join(output_directory, "Runtime_log.txt"))

        description =  "<mapDamage: %i file(s) -> '%s'>" % (len(input_files), output_directory)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([cmd_cat, cmd_map]),
                             description  = description,
                             dependencies = dependencies)


class FilterUniqueBAMNode(CommandNode):
    def __init__(self, config, input_bams, output_bam, dependencies = ()):
        merge_jar  = os.path.join(config.jar_root, "MergeSamFiles.jar")
        merge_call = ["java", "-jar", merge_jar, 
                      "TMP_DIR=%s" % config.temp_root, 
                      "SO=coordinate",
                      "QUIET=true",
                      "COMPRESSION_LEVEL=0",
                      "OUTPUT=/dev/stdout"]
        merge_files = {"OUT_STDOUT" : AtomicCmd.PIPE}
        for (ii, filename) in enumerate(input_bams, start = 1):
            merge_call.append("INPUT=%%(IN_FILE_%i)s" % ii)
            merge_files["IN_FILE_%i" % ii] = filename


        merge = AtomicCmd(merge_call, **merge_files)
        filteruniq = AtomicCmd(["FilterUniqueBAM", "--PIPE", "--library"],
                               IN_STDIN   = merge,
                               OUT_STDOUT = output_bam)

        command     = ParallelCmds([merge, filteruniq])
        description =  "<FilterUniqueBAM: %s>" % (self._desc_files(input_bams),)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


class ValidateBAMFile(MetaNode):
    def __init__(self, config, node, log_file = None, dependencies = None):
        filenames = self._get_input_file(node)
        assert len(filenames) == 1, (filenames, node)
        input_file = filenames.pop()

        subnodes = [node]
        if isinstance(node, MetaNode):
            subnodes = list(node.subnodes)

        validation_params = ValidateBAMNode.customize(config       = config,
                                                      input_bam    = input_file,
                                                      output_log   = log_file,
                                                      dependencies = subnodes)
        # Ignored since we filter out misses and low-quality hits during mapping, which
        # leads to a large proportion of missing mates for PE reads.
        validation_params.command.push_parameter("IGNORE", "MATE_NOT_FOUND", sep = "=")
        # Ignored due to high rate of false positives for lanes with few hits, where
        # high-quality reads may case ValidateSamFile to mis-identify the qualities
        validation_params.command.push_parameter("IGNORE", "INVALID_QUALITY_FORMAT", sep = "=")

        subnodes.append(validation_params.build_node())

        description = "<w/Validation: " + str(node)[1:]
        MetaNode.__init__(self,
                          description  = description,
                          subnodes     = subnodes,
                          dependencies = dependencies or node.dependencies)


    def _get_input_file(cls, node):
        filenames = set()
        for filename in node.output_files:
            if filename.lower().endswith(".bai"):
                filenames.add(swap_ext(filename, ".bam"))
            elif filename.lower().endswith(".bam"):
                filenames.add(filename)

        if not filenames and node.subnodes:
            for subnode in node.subnodes:
                filenames.update(cls._get_input_file(subnode))

        return filenames


class CleanupBAMNode(CommandNode):
    def __init__(self, config, reference, input_bam, output_bam, min_mapq, tags, dependencies = ()):
        flt = AtomicCmd(["samtools", "view", "-bu", "-F0x4", "-q%i" % min_mapq, "%(IN_BAM)s"],
                        IN_BAM  = input_bam,
                        OUT_STDOUT = AtomicCmd.PIPE)

        jar_file = os.path.join(config.jar_root, "AddOrReplaceReadGroups.jar")
        params = AtomicJavaParams(config, jar_file)
        params.set_parameter("INPUT", "/dev/stdin", sep = "=")
        params.set_parameter("OUTPUT", "/dev/stdout", sep = "=")
        params.set_parameter("QUIET", "true", sep = "=")
        params.set_parameter("COMPRESSION_LEVEL", "0", sep = "=")

        for (tag, value) in sorted(tags.iteritems()):
            params.set_parameter(tag, value, sep = "=")

        params.set_paths(IN_STDIN   = flt,
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


def IndexedFilterUniqueBAMNode(output_bam, **kwargs):
    node = FilterUniqueBAMNode(output_bam = output_bam,
                               **kwargs)

    return BAMIndexNode(infile       = output_bam,
                        dependencies = node)
