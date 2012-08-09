import os

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds
from pypeline.tools import sam_to_bam



class SE_BWANode(CommandNode):
    def __init__(self, input_file, output_file, prefix, read_group, threads = 1, dependencies = ()):
        aln   = AtomicCmd(["bwa", "aln", 
                           "-l", 2 ** 20, 
                           "-t", threads,
                           prefix, "%(IN_FILE)s"],
                          IN_FILE = input_file,
                          OUT_STDOUT = AtomicCmd.PIPE)
        
        samse = AtomicCmd(["bwa", "samse", "-r", read_group, prefix, "-", "%(IN_FILE)s"],
                          IN_STDIN = aln,
                          IN_FILE  = input_file,
                          OUT_STDOUT = AtomicCmd.PIPE)

        flt = sam_to_bam.build_atomiccmd(AtomicCmd,
                                         min_quality   = 25,
                                         exclude_flags = 0x4,
                                         stdin         = samse,
                                         output_file   = output_file)

        description =  "<SE_BWA (%i threads): '%s' -> '%s'>" % (threads, input_file, output_file)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln, samse, flt]),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)


class PE_BWANode(CommandNode):
    def __init__(self, input_file_1, input_file_2, output_file, prefix, read_group, threads = 1, dependencies = ()):
        aln_call = ["bwa", "aln", 
                    "-l", 2 ** 20, 
                    "-t", max(1, threads / 2),
                    prefix, "%(IN_FILE)s", "-f", "%(TEMP_OUT_SAI)s"]
        aln_1 = AtomicCmd(aln_call,
                          IN_FILE  = input_file_1,
                          TEMP_OUT_SAI = "pair_1.sai")

        aln_2 = AtomicCmd(aln_call,
                          IN_FILE  = input_file_2,
                          TEMP_OUT_SAI = "pair_2.sai")
        
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
                                         min_quality   = 25,
                                         exclude_flags = 0x4,
                                         stdin         = samse,
                                         output_file   = output_file)

        description =  "<PE_BWA (%i threads): '%s' -> '%s'>" % (threads, input_file_1, output_file)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln_1, aln_2, samse, flt]),
                             description  = description,
                             threads      = threads,
                             dependencies = dependencies)


    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "pair_1.sai"))
        os.mkfifo(os.path.join(temp, "pair_2.sai"))
