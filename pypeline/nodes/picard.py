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
from pypeline.atomicparams import *
from pypeline.common.fileutils import swap_ext
from pypeline.common.utilities import safe_coerce_to_tuple


class ValidateBAMNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bam, output_log = None, dependencies = ()):
        jar_file = os.path.join(config.picard_root, "ValidateSamFile.jar")
        params = AtomicJavaParams(config, jar_file)

        params.set_parameter("I", "%(IN_BAM)s", sep = "=")
        params.set_paths(IN_BAM     = input_bam,
                         OUT_STDOUT = output_log or (input_bam + ".log"))
       
        return {"command" : params}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        CommandNode.__init__(self, 
                             command      = parameters.command.create_cmd(),
                             description  = "<Validate BAM: '%s'>" % (parameters.input_bam,),
                             dependencies = parameters.dependencies)


class MarkDuplicatesNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bams, output_bam, output_metrics = None, dependencies = ()):
        jar_file = os.path.join(config.picard_root, "MarkDuplicates.jar")
        params = AtomicJavaParams(config, jar_file)

        # Create .bai by default, since it is required by a lot of other programs
        params.set_parameter("CREATE_INDEX", "True", sep = "=", fixed = False)
        # Remove duplicates from output by default to save disk-space
        params.set_parameter("REMOVE_DUPLICATES", "True", sep = "=", fixed = False)

        params.set_option("OUTPUT", "%(OUT_BAM)s", sep = "=")
        params.set_option("METRICS_FILE", "%(OUT_METRICS)s", sep = "=")

        for (index, filename) in enumerate(safe_coerce_to_tuple(input_files)):
            params.push_parameter("I", "%%(IN_BAM_%02i)s" % index)
            params.set_paths("IN_BAM_%02i" % index, filename)

        params.set_paths(OUT_BAM     = output_file,
                         OUT_BAI     = swap_ext(output_file, ".bai"),
                         OUT_METRICS = metrics_file or (output_file + ".metrics"))

        return {"command" : params}
        

    def __init__(self, parameters):
        description =  "<MarkDuplicates: %s>" % (self._desc_files(parameters.input_bams),)
        CommandNode.__init__(self, 
                             command      = parameters.command.create_cmd(),
                             description  = description,
                             dependencies = parameters.dependencies)


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


if __name__ == '__main__':
    class Config:
        temp_root = "temp/tmp"
        picard_root = "/home/mischu/archive/research/tools/bamPipeline/picard-tools-1.69"

    params = ValidateBAMNode.customize(config     = Config, 
                                       input_bam  = "/home/mischu/archive/research/tools/seqStats/test/data/exampleBAM.bam",
                                       output_log = "temp/validated.log")
    params.command.push_parameter("IGNORE=MATE_NOT_FOUND")
    params.command.push_parameter("IGNORE=MISSING_TAG_NM")
    node = params.build_node()

    node.run(Config)
