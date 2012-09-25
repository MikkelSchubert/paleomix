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
from pypeline.common.utilities import safe_coerce_to_tuple


class ValidateBAMNode(CommandNode):
    def __init__(self, config, bamfile, output_file = None, ignore = (), dependencies = ()):
        jar  = os.path.join(config.picard_root, "ValidateSamFile.jar")
        call = ["java", "-jar", jar,
                "TMP_DIR=%s" % config.temp_root, "INPUT=%(IN_BAM)s"]

        for error in ignore:
            call.append("IGNORE=%s" % (error,))

        if not output_file:
            output_file = bamfile + ".validated"

        command = AtomicCmd(call,
                            IN_BAM     = bamfile,
                            OUT_STDOUT = output_file)

        CommandNode.__init__(self, 
                             command      = command,
                             description  = "<Validate BAM: '%s'>" % (bamfile,),
                             dependencies = dependencies)


class MarkDuplicatesNode(CommandNode):
    def __init__(self, config, input_files, output_file, metrics_file = None, keep_duplicates = False, dependencies = ()):
        if not metrics_file:
            metrics_file = output_file + ".metrics"

        jar  = os.path.join(config.picard_root, "MarkDuplicates.jar")
        call = ["java", "-jar", jar, 
                "TMP_DIR=%s" % config.temp_root, 
                "REMOVE_DUPLICATES=%s" % str(not keep_duplicates).lower(),
                "CREATE_INDEX=True",
                # ASSUME_SORTED is required to allow use of 'samtools sort', which
                # does not add a line to the header specifying sorting. In the case
                # of unsorted files the command will abort.
                "ASSUME_SORTED=True", 
                "OUTPUT=%(OUT_BAM)s",
                "METRICS_FILE=%(OUT_METRICS)s"]

        args = {"OUT_BAM"     : output_file,
                "OUT_BAI"     : swap_ext(output_file, ".bai"),
                "OUT_METRICS" : metrics_file}       

        input_files = safe_coerce_to_tuple(input_files)
        for (index, input_file) in enumerate(input_files):
            call.append("INPUT=%%(IN_BAM_%i)s" % index)
            args["IN_BAM_%i" % index] = input_file

        description =  "<MarkDuplicates: %s>" % (self._desc_files(input_files),)
        CommandNode.__init__(self, 
                             command      = AtomicCmd(call, **args),
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
