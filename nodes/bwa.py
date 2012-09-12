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

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds
from pypeline.tools import sam_to_bam


class BWAIndex(CommandNode):
    def __init__(self, input_file, prefix = None, algorithm = "is", dependencies = ()):
        prefix = prefix if prefix else input_file
        command = AtomicCmd(["bwa", "index", 
                             "-a", algorthm,
                             "-p", prefix]
                            IN_FILE = input_file,
                            **_prefix_files(prefix, iotype = "OUT"))
        
        description =  "<BWA Index '%s' -> '%s'>" % (threads, input_file, prefix)
        CommandNode.__init__(self, 
                             command      = command
                             description  = description,
                             dependencies = dependencies)


class SE_BWANode(CommandNode):
    def __init__(self, input_file, output_file, prefix, read_group, min_quality = 0, threads = 1, dependencies = ()):
        aln   = AtomicCmd(["bwa", "aln", 
                           "-l", 2 ** 20, 
                           "-t", threads,
                           prefix, "%(IN_FILE)s"],
                          IN_FILE = input_file,
                          OUT_STDOUT = AtomicCmd.PIPE,
                          **_prefix_files(prefix))
        
        samse = AtomicCmd(["bwa", "samse", "-r", read_group, prefix, "-", "%(IN_FILE)s"],
                          IN_STDIN = aln,
                          IN_FILE  = input_file,
                          OUT_STDOUT = AtomicCmd.PIPE)

        flt = sam_to_bam.build_atomiccmd(AtomicCmd,
                                         min_quality   = min_quality,
                                         exclude_flags = 0x4,
                                         stdin         = samse,
                                         output_file   = AtomicCmd.PIPE)

        sort = AtomicCmd(["samtools", "sort", "-", "%(TEMP_OUT_BAM)s"],
                         IN_STDIN     = flt,
                         OUT_BAM      = output_file,
                         # Prefix used by 'samtools sort', with .bam added to final file
                         TEMP_OUT_BAM = os.path.splitext(output_file)[0])


        description =  "<SE_BWA (%i threads): '%s'>" % (threads, input_file)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln, samse, flt, sort]),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)


class PE_BWANode(CommandNode):
    def __init__(self, input_file_1, input_file_2, output_file, prefix, read_group, min_quality = 0, threads = 1, dependencies = ()):
        aln_call = ["bwa", "aln", 
                    "-l", 2 ** 20, 
                    "-t", max(1, threads / 2),
                    prefix, "%(IN_FILE)s", "-f", "%(TEMP_OUT_SAI)s"]
        aln_1 = AtomicCmd(aln_call,
                          IN_FILE  = input_file_1,
                          TEMP_OUT_SAI = "pair_1.sai",
                          **_prefix_files(prefix))

        aln_2 = AtomicCmd(aln_call,
                          IN_FILE  = input_file_2,
                          TEMP_OUT_SAI = "pair_2.sai",
                          **_prefix_files(prefix))
        
        samse = AtomicCmd(["bwa", "sampe", prefix, 
                           "-P", "-r", read_group,
                           "%(TEMP_IN_SAI_1)s", "%(TEMP_IN_SAI_2)s", 
                           "%(IN_FILE_1)s",     "%(IN_FILE_2)s"],
                          IN_FILE_1     = input_file_1,
                          IN_FILE_2     = input_file_2,
                          TEMP_IN_SAI_1 = "pair_1.sai",
                          TEMP_IN_SAI_2 = "pair_2.sai",
                          OUT_STDOUT    = AtomicCmd.PIPE)

        flt = sam_to_bam.build_atomiccmd(AtomicCmd,
                                         min_quality   = min_quality,
                                         exclude_flags = 0x4,
                                         stdin         = samse,
                                         output_file   = AtomicCmd.PIPE)

        sort = AtomicCmd(["samtools", "sort", "-", "%(TEMP_OUT_BAM)s"],
                         IN_STDIN     = flt,
                         OUT_BAM      = output_file,
                         # Prefix used by 'samtools sort', with .bam added to final file
                         TEMP_OUT_BAM = os.path.splitext(output_file)[0])


        description =  "<PE_BWA (%i threads): '%s'>" % (threads, input_file_1)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln_1, aln_2, samse, flt, sort]),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)


    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "pair_1.sai"))
        os.mkfifo(os.path.join(temp, "pair_2.sai"))


class BWASWNode(CommandNode):
    def __init__(self, input_file_1, output_file, prefix, input_file_2 = None, min_quality = 0, threads = 1, dependencies = ()):
        command = ["bwa", "bwasw", "-t", threads, prefix, "%(IN_FILE_1)s"]
        files   = {"IN_FILE_1"  : input_file_1,
                   "OUT_STDOUT" : AtomicCmd.PIPE}
        if input_file_2:
            command.append("%(IN_FILE_2)s")
            files["IN_FILE_2"] = input_file_2
        files.update(_prefix_files(prefix))

        aln = AtomicCmd(command, **files)
        flt = sam_to_bam.build_atomiccmd(AtomicCmd,
                                         min_quality   = min_quality,
                                         exclude_flags = 0x4,
                                         stdin         = aln,
                                         output_file   = output_file)

        description =  "<PE_BWASW (%i threads): '%s' -> '%s'>" % (threads, input_file_1, output_file)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln, flt]),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)



def _prefix_files(prefix, iotype = "IN"):
    files = {}
    for postfix in ("amb", "ann", "bwt", "pac", "rbwt", "rpac", "rsa", "sa"):
        files["%s_PREFIX_%s" % (iotype, postfix.upper())] = prefix + "." + postfix
    return files
