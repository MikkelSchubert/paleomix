import os

import node
import fileutils

from atomiccmd import AtomicCmd



class _IndelTrainerNode(node.Node):
    def __init__(self, config, destination, reference, infile, outfile, dependencies = ()):
        jarfile = os.path.join(config.gatk_root, "GenomeAnalysisTK.jar")
        outfile += ".intervals"
        call  = ["java", "-jar", "%(IN_JAR)s",
                 "-T", "RealignerTargetCreator", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-o", "%(OUT_INTERVALS)s"]

        command = AtomicCmd(destination, 
                            call, 
                            IN_JAR = jarfile,
                            IN_REFERENCE = reference,
                            IN_BAMFILE = infile,
                            OUT_INTERVALS = outfile,
                            stdout = outfile + ".log",
                            stderr = outfile + ".log")
        
        description = "<Train Indel Realigner: '%s' -> '%s'>" \
            % (infile, os.path.join(destination, outfile))

        node.Node.__init__(self, 
                           description = description,
                           command = command,
                           dependencies = dependencies)


class _IndelRealignerNode(node.Node):
    def __init__(self, config, destination, reference, infile, outfile, dependencies = ()):
        jarfile = os.path.join(config.gatk_root, "GenomeAnalysisTK.jar")
        intervalsfile = os.path.join(destination, outfile + ".intervals")
        call  = ["java", "-jar", "%(IN_JAR)s",
                 "-T", "IndelRealigner", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-targetIntervals", "%(IN_INTERVALS)s",
                 "-o", "%(OUT_BAMFILE)s"]

        command = AtomicCmd(destination,
                            call, 
                            IN_JAR = jarfile,
                            IN_REFERENCE = reference,
                            IN_BAMFILE   = infile,
                            IN_INTERVALS = intervalsfile,
                            OUT_BAMFILE = outfile,
                            OUT_INDEX  = fileutils.swap_ext(outfile, ".bai"),
                            stdout = outfile + ".log",
                            stderr = outfile + ".log")

        description = "<Indel Realigner: '%s' -> '%s'>" \
            % (infile, os.path.join(destination, outfile))

        node.Node.__init__(self, 
                           description = description,
                           command = command,
                           dependencies = dependencies)



def IndelRealignerNode(config, destination, reference, infile, outfile, dependencies = ()):
    trainer = _IndelTrainerNode(config = config,
                                destination = destination, 
                                reference = reference, 
                                infile = infile,
                                outfile = outfile,
                                dependencies = dependencies)
    aligner = _IndelRealignerNode(config = config,
                                  destination = destination, 
                                  reference = reference, 
                                  infile = infile, 
                                  outfile = outfile,
                                  dependencies = [trainer])
    return aligner
