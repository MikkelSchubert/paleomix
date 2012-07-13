import os

import pypeline.fileutils as fileutils
from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd import AtomicCmd



class _IndelTrainerNode(CommandNode):
    def __init__(self, config, reference, infile, outfile, dependencies = ()):
        jarfile = os.path.join(config.gatk_root, "GenomeAnalysisTK.jar")
        outfile += ".intervals"
        call  = ["java", "-jar", jarfile,
                 "-T", "RealignerTargetCreator", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-o", "%(OUT_INTERVALS)s"]

        command = AtomicCmd(call, 
                            IN_REFERENCE  = reference,
                            IN_BAMFILE    = infile,
                            OUT_INTERVALS = outfile,
                            OUT_STDOUT    = outfile + ".log")
        
        description = "<Train Indel Realigner: '%s' -> '%s'>" \
            % (infile, outfile)

        CommandNode.__init__(self, 
                             description = description,
                             command = command,
                             dependencies = dependencies)



class _IndelRealignerNode(CommandNode):
    def __init__(self, config, reference, infile, outfile, dependencies = ()):
        jarfile = os.path.join(config.gatk_root, "GenomeAnalysisTK.jar")
        intervalsfile = outfile + ".intervals"
        call  = ["java", "-jar", jarfile,
                 "-T", "IndelRealigner", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-targetIntervals", "%(IN_INTERVALS)s",
                 "-o", "%(OUT_BAMFILE)s"]
        
        command = AtomicCmd(call, 
                            IN_REFERENCE = reference,
                            IN_BAMFILE   = infile,
                            IN_INTERVALS = intervalsfile,
                            OUT_BAMFILE = outfile,
                            OUT_INDEX  = fileutils.swap_ext(outfile, ".bai"),
                            OUT_STDOUT = outfile + ".log")

        description = "<Indel Realign: '%s' -> '%s'>" \
            % (infile, outfile)

        CommandNode.__init__(self, 
                             description = description,
                             command = command,
                             dependencies = dependencies)



class IndelRealignerNode(MetaNode):
    def __init__(self, config, reference, infile, outfile, dependencies = ()):
        trainer = _IndelTrainerNode(config = config,
                                    reference = reference, 
                                    infile = infile,
                                    outfile = outfile,
                                    dependencies = dependencies)
        aligner = _IndelRealignerNode(config = config,
                                      reference = reference, 
                                      infile = infile, 
                                      outfile = outfile,
                                      dependencies = trainer)
        
        MetaNode.__init__(self, 
                          description  = "<GATK Indel Realigner: '%s' -> '%s'>" \
                              % (infile, outfile),
                          subnodes     = [trainer, aligner],
                          dependencies = dependencies)
