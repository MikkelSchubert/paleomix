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


class BWAIndexNode(CommandNode):
    def __init__(self, input_file, prefix = None, algorithm = "is", dependencies = ()):
        prefix = prefix if prefix else input_file
        command = AtomicCmd(["bwa", "index", 
                             "-a", algorithm,
                             "-p", "%(TEMP_OUT_PREFIX)s",
                             "%(IN_FILE)s"],
                            IN_FILE = input_file,
                            TEMP_OUT_PREFIX = os.path.basename(prefix),
                            **_prefix_files(prefix, iotype = "OUT"))
        
        description =  "<BWA Index '%s' -> '%s.*'>" % (input_file, prefix)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


class SE_BWANode(CommandNode):
    def __init__(self, input_file, output_file, reference, prefix, read_group, min_quality = 0, threads = 1, dependencies = ()):
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

        cmds = _process_output(samse, output_file, reference, min_quality)

        description =  "<SE_BWA (%i threads): '%s'>" % (threads, input_file)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln, samse] + cmds),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)


class PE_BWANode(CommandNode):
    def __init__(self, input_file_1, input_file_2, output_file, reference, prefix, read_group, min_quality = 0, threads = 1, dependencies = ()):
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

        cmds = _process_output(samse, output_file, reference, min_quality)

        description =  "<PE_BWA (%i threads): '%s'>" % (threads, input_file_1)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln_1, aln_2, samse] + cmds),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)


    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "pair_1.sai"))
        os.mkfifo(os.path.join(temp, "pair_2.sai"))


class BWASWNode(CommandNode):
    def __init__(self, input_file_1, output_file, reference, prefix, input_file_2 = None, min_quality = 0, threads = 1, parameters = [], dependencies = ()):
        command = ["bwa", "bwasw", "-t", threads, prefix, "%(IN_FILE_1)s"] + parameters
        files   = {"IN_FILE_1"  : input_file_1,
                   "OUT_STDOUT" : AtomicCmd.PIPE}
        if input_file_2:
            command.append("%(IN_FILE_2)s")
            files["IN_FILE_2"] = input_file_2
        files.update(_prefix_files(prefix))

        aln = AtomicCmd(command, **files)
        cmds = _process_output(samse, output_file, reference, min_quality)


        if input_file_2:
            description =  "<PE_BWASW (%i threads): '%s', '%s' -> '%s'>" % (threads, input_file_1, input_file_2, output_file)
        else:
            description =  "<BWASW (%i threads): '%s' -> '%s'>" % (threads, input_file_1, output_file)

        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln] + cmds),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)


def _process_output(stdin, output_file, reference, min_quality):
    convert = AtomicCmd(["safeSAM2BAM", "--flag-as-sorted"],
                        IN_STDIN   = stdin,
                        OUT_STDOUT = AtomicCmd.PIPE)

    flt = AtomicCmd(["samtools", "view", "-bu" "-F0x4", "-q%i" % min_quality, "-"],
                    IN_STDIN  = convert,
                    OUT_STDOUT = AtomicCmd.PIPE)

    sort = AtomicCmd(["samtools", "sort", "-o", "-", "%(TEMP_OUT_BAM)s"],
                     IN_STDIN     = flt,
                     OUT_STDOUT   = AtomicCmd.PIPE,
                     TEMP_OUT_BAM = "sorted")

    calmd = AtomicCmd(["samtools", "calmd", "-b", "-", "%(IN_REF)s"],
                      IN_REF   = reference,
                      IN_STDIN = sort,
                      OUT_STDOUT = output_file)

    return [convert, flt, sort, calmd]



def _prefix_files(prefix, iotype = "IN"):
    files = {}
    for postfix in ("amb", "ann", "bwt", "pac", "rbwt", "rpac", "rsa", "sa"):
        files["%s_PREFIX_%s" % (iotype, postfix.upper())] = prefix + "." + postfix
    return files
