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
from pypeline.common.fileutils import swap_ext


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
                "REMOVE_DUPLICATES=%s" % str(not keep_duplicates).lower(),
                "CREATE_INDEX=True",
                "INPUT=%(IN_BAM)s",
                "OUTPUT=%(OUT_BAM)s",
                "METRICS_FILE=%(OUT_METRICS)s"]

        command = AtomicCmd(call,
                            IN_BAM      = input_file,
                            OUT_BAM     = output_file,
                            OUT_BAI     = swap_ext(output_file, ".bai"),
                            OUT_METRICS = output_file + ".metrics")

        description =  "<MarkDuplicates: '%s' -> '%s'>" % (input_file, output_file)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


class MergeSamFilesNode(CommandNode):
    def __init__(self, config, input_files, output_file, create_index = True, dependencies = ()):
        jar  = os.path.join(config.picard_root, "MergeSamFiles.jar")
        call = ["java", "-jar", jar, 
                "TMP_DIR=%s" % config.temp_root, 
                "SO=coordinate",
                "OUTPUT=%(OUT_BAM)s"]
        files = {"OUT_BAM" : output_file}

        if create_index:
            call.append("CREATE_INDEX=True")
            files["OUT_BAI"] = swap_ext(output_file, ".bai")


        for (ii, filename) in enumerate(input_files, start = 1):
            call.append("INPUT=%%(IN_FILE_%i)s" % ii)
            files["IN_FILE_%i" % ii] = filename
        
        command = AtomicCmd(call, **files)

        description =  "<Merge BAMs: %i file(s) -> '%s'>" % (len(input_files), output_file)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = dependencies)
