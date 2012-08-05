import os

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd


class ValidateBAMNode(CommandNode):
    def __init__(self, config, bamfile, ignore = (), dependencies = ()):
        jar  = os.path.join(config.picard_root, "ValidateSamFile.jar")
        call = ["java", "-jar", jar,
                "TMP_DIR=%s" % config.temp_root, "INPUT=%(IN_BAM)s"]
        for error in ignore:
            call.append("IGNORE=%s" % (error,))

        command = AtomicCmd(call,
                            IN_BAM     = bamfile,
                            OUT_STDOUT = bamfile + ".validated")

        CommandNode.__init__(self, 
                             command      = command,
                             description  = "<Validate BAM: '%s'>" % (bamfile,),
                             dependencies = dependencies)


class MarkDuplicatesNode(CommandNode):
    def __init__(self, config, input_file, output_file, keep_duplicates = False, dependencies = ()):
        jar  = os.path.join(config.picard_root, "MarkDuplicates.jar")
        call = ["java", "-jar", jar, 
                "TMP_DIR=%s" % config.temp_root, 
                "REMOVE_DUPLICATES=%s" % str(keep_duplicates).lower(),
                "INPUT=%(IN_BAM)s",
                "OUTPUT=%(OUT_BAM)s",
                "METRICS_FILE=%(OUT_METRICS)s"]

        command = AtomicCmd(call,
                            IN_BAM      = input_file,
                            OUT_BAM     = output_file,
                            OUT_METRICS = output_file + ".metrics")

        description =  "<MarkDuplicates: '%s' -> '%s'>" % (input_file, output_file)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


class MergeSamFilesNode(CommandNode):
    def __init__(self, config, input_files, output_file, dependencies = ()):
        jar  = os.path.join(config.picard_root, "MergeSamFiles.jar")
        call = ["java", "-jar", jar, 
                "TMP_DIR=%s" % config.temp_root, 
                "SO=coordinate",
                "OUTPUT=%(OUT_BAM)s"]
        files = {"OUT_BAM" : output_file}

        for (ii, filename) in enumerate(input_files, start = 1):
            call.append("INPUT=%%(IN_FILE_%i)s" % ii)
            files["IN_FILE_%i" % ii] = filename
        
        command = AtomicCmd(call, **files)

        description =  "<Merge BAMs: %i file(s) -> '%s'>" % (len(input_files), output_file)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)
