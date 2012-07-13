import os

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd


class ValidateBAMNode(CommandNode):
    def __init__(self, config, bamfile, ignore = (), dependencies = ()):
        call  = ["java", "-jar", "%(IN_JAR)s", "TMP_DIR=%s" % config.temp_root, "I=%(IN_BAM)s"]
        for error in ignore:
            call.append("IGNORE=%s" % (error,))

        jarfile = os.path.join(config.picard_root, "ValidateSamFile.jar")
        logfile = bamfile + ".validated"

        command = AtomicCmd(call,
                            IN_JAR = jarfile,
                            IN_BAM = bamfile,
                            OUT_STDOUT = logfile,
                            OUT_STDERR = logfile)

        description =  "<Validate BAM: '%s' -> '%s'>" % (bamfile, logfile)

        CommandNode.__init__(self, 
                             description = description,
                             command = command,
                             dependencies = dependencies)
