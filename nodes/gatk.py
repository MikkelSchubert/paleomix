import os

import node
import fileutils

from atomiccmd import AtomicCmd



class _IndelTrainerNode(node.Node):
    def __init__(self, config, destination, reference, inFile, outFile, dependencies = []):
        jarFile = os.path.join(config.gatk_root, "GenomeAnalysisTK.jar")
        outFile += ".intervals"
        call  = ["java", "-jar", "%(IN_JAR)s",
                 "-T", "RealignerTargetCreator", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-o", "%(OUT_INTERVALS)s"]

        command = AtomicCmd(destination, 
                            call, 
                            IN_JAR = jarFile,
                            IN_REFERENCE = reference,
                            IN_BAMFILE = inFile,
                            OUT_INTERVALS = outFile,
                            stdout = outFile + ".log",
                            stderr = outFile + ".log")
        
        description = "<Train Indel Realigner: '%s' -> '%s'>" \
            % (inFile, os.path.join(destination, outFile))

        node.Node.__init__(self, 
                           description = description,
                           command = command,
                           dependencies = dependencies)


class _IndelRealignerNode(node.Node):
    def __init__(self, config, destination, reference, inFile, outFile, dependencies = []):
        jarFile = os.path.join(config.gatk_root, "GenomeAnalysisTK.jar")
        intervalsFile = os.path.join(destination, outFile + ".intervals")
        call  = ["java", "-jar", "%(IN_JAR)s",
                 "-T", "IndelRealigner", 
                 "-R", "%(IN_REFERENCE)s",
                 "-I", "%(IN_BAMFILE)s",
                 "-targetIntervals", "%(IN_INTERVALS)s",
                 "-o", "%(OUT_BAMFILE)s"]

        command = AtomicCmd(destination,
                            call, 
                            IN_JAR = jarFile,
                            IN_REFERENCE = reference,
                            IN_BAMFILE   = inFile,
                            IN_INTERVALS = intervalsFile,
                            OUT_BAMFILE = outFile,
                            OUT_INDEX  = fileutils.swap_ext(outFile, ".bai"),
                            stdout = outFile + ".log",
                            stderr = outFile + ".log")

        description = "<Indel Realigner: '%s' -> '%s'>" \
            % (inFile, os.path.join(destination, outFile))

        node.Node.__init__(self, 
                           description = description,
                           command = command,
                           dependencies = dependencies)



def add_indel_realigner(config, destination, reference, inFile, outFile, dependencies = []):
    trainer = _IndelTrainerNode(config = config,
                                destination = destination, 
                                reference = reference, 
                                inFile = inFile,
                                outFile = outFile,
                                dependencies = dependencies)
    aligner = _IndelRealignerNode(config = config,
                                  destination = destination, 
                                  reference = reference, 
                                  inFile = inFile, 
                                  outFile = outFile,
                                  dependencies = [trainer])

    return aligner
