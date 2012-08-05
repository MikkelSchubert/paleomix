import os

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds


class SE_BWANode(CommandNode):
    def __init__(self, input_file, output_file, prefix, read_group, dependencies = ()):
        aln   = AtomicCmd(["bwa", "aln", "-l", 2 ** 20, prefix, "%(IN_FILE)s"],
                          IN_FILE = input_file,
                          OUT_STDOUT = AtomicCmd.PIPE)
        
        samse = AtomicCmd(["bwa", "samse", "-r", read_group, prefix, "-", "%(IN_FILE)s"],
                          IN_STDIN = aln,
                          IN_FILE  = input_file,
                          OUT_STDOUT = AtomicCmd.PIPE)

        flt = AtomicCmd(["samtools", "view", "-bSh", "-q25", "-F0x4", "-"],
                        IN_STDIN   = samse,
                        OUT_STDOUT = output_file)

        description =  "<SE_BWA: '%s' -> '%s'>" % (input_file, output_file)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln, samse, flt]),
                             description  = description,
                             dependencies = dependencies)


class PE_BWANode(CommandNode):
    def __init__(self, input_file_1, input_file_2, output_file, prefix, read_group, dependencies = ()):
        aln_call = ["bwa", "aln", "-l", 2 ** 20, prefix, "%(IN_FILE)s", "-f", "%(TEMP_OUT_SAI)s"]
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

        flt = AtomicCmd(["samtools", "view", "-bSh", "-q25", "-F0x4", "-"],
                        IN_STDIN   = samse,
                        OUT_STDOUT = output_file)

        description =  "<PE_BWA: '%s' -> '%s'>" % (input_file_1, output_file)
        CommandNode.__init__(self, 
                             command      = ParallelCmds([aln_1, aln_2, samse, flt]),
                             description  = description,
                             dependencies = dependencies)


    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "pair_1.sai"))
        os.mkfifo(os.path.join(temp, "pair_2.sai"))
