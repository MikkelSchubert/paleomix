import os

import node
from atomiccmd import AtomicCmd


class ValidateBAMNode(node.SimpleNode):
    def __init__(self, config, bamfile, ignore = (), dependencies = ()):
        call  = ["java", "-jar", "%(IN_JAR)s", "TMP_DIR=%s" % config.temp_root, "I=%(IN_BAM)s"]
        for error in ignore:
            call.append("IGNORE=%s" % (error,))

        jarfile = os.path.join(config.jar_root, "ValidateSamFile.jar")
        logfile = os.path.split(bamfile)[-1] + ".validated"
        destination = os.path.split(bamfile)[0]

        command = AtomicCmd(destination,
                              call,
                              IN_JAR = jarfile,
                              IN_BAM = bamfile,
                              stdout = logfile,
                              stderr = logfile)

        description =  "<Validate BAM: '%s' -> '%s'>" \
            % (bamfile, os.path.join(destination, logfile))

        node.SimpleNode.__init__(self, 
                                 description = description,
                                 command = command,
                                 dependencies = dependencies)
