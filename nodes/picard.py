import os

import node
import fileutils
from atomiccmd import AtomicCmd


class ValidateBAMNode(node.Node):
    def __init__(self, config, bamFile, ignore = [], dependencies = []):
        call  = ["java", "-jar", "%(IN_JAR)s", "I=%(IN_BAM)s"]
        for error in ignore:
            call.append("IGNORE=%s" % (error,))

        jarFile = os.path.join(config.jar_root, "ValidateSamFile.jar")
        logFile = os.path.split(bamFile)[-1] + ".validated"
        destination = os.path.split(bamFile)[0]

        command = AtomicCmd(destination,
                              call,
                              IN_JAR = jarFile,
                              IN_BAM = bamFile,
                              stdout = logFile,
                              stderr = logFile)

        description =  "<Validate BAM: '%s' -> '%s'>" \
            % (bamFile, os.path.join(destination, logFile))

        node.Node.__init__(self, 
                           description = description,
                           command = command,
                           dependencies = dependencies)
