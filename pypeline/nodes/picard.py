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
        jar_file = os.path.join(config.jar_root, "ValidateSamFile.jar")
        params = AtomicJavaParams(config, jar_file)

        params.set_parameter("I", "%(IN_BAM)s", sep = "=")
        params.set_paths(IN_BAM     = input_bam,
                         OUT_STDOUT = output_log or swap_ext(input_bam, ".validated"))

        return {"command" : params}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<Validate BAM: '%s'>" % (parameters.input_bam,),
                             dependencies = parameters.dependencies)


class BuildSequenceDictNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, config, reference, dependencies = ()):
        jar_file = os.path.join(config.jar_root, "CreateSequenceDictionary.jar")
        params = AtomicJavaParams(config, jar_file)

        params.set_parameter("R", "%(IN_REF)s", sep = "=")
        params.set_parameter("O", "%(OUT_DICT)s", sep = "=")
        params.set_paths(IN_REF     = reference,
                         OUT_DICT   = swap_ext(reference, ".dict"))

        return {"command" : params}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = "<SequenceDictionary: '%s'>" % (parameters.reference,),
                             dependencies = parameters.dependencies)



class MarkDuplicatesNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bams, output_bam, output_metrics = None, dependencies = ()):
        jar_file = os.path.join(config.jar_root, "MarkDuplicates.jar")
        params = AtomicJavaParams(config, jar_file)

        # Create .bai index, since it is required by a lot of other programs
        params.set_parameter("CREATE_INDEX", "True", sep = "=")

        params.set_parameter("OUTPUT", "%(OUT_BAM)s", sep = "=")
        params.set_parameter("METRICS_FILE", "%(OUT_METRICS)s", sep = "=")

        input_bams = safe_coerce_to_tuple(input_bams)
        for (index, filename) in enumerate(input_bams):
            params.push_parameter("I", "%%(IN_BAM_%02i)s" % index, sep = "=")
            params.set_paths("IN_BAM_%02i" % index, filename)

        # Remove duplicates from output by default to save disk-space
        params.set_parameter("REMOVE_DUPLICATES", "True", sep = "=", fixed = False)

        params.set_paths(OUT_BAM     = output_bam,
                         OUT_BAI     = swap_ext(output_bam, ".bai"),
                         OUT_METRICS = output_metrics or swap_ext(output_bam, ".metrics"))

        return {"command" : params}
        

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        description =  "<MarkDuplicates: %s>" % (self._desc_files(parameters.input_bams),)
        CommandNode.__init__(self, 
                             command      = parameters.command.finalize(),
                             description  = description,
                             dependencies = parameters.dependencies)




class MergeSamFilesNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bams, output_bam, dependencies = ()):
        jar_file = os.path.join(config.jar_root, "MergeSamFiles.jar")
        params = AtomicJavaParams(config, jar_file)
        
        params.set_parameter("OUTPUT", "%(OUT_BAM)s", sep = "=")
        params.set_parameter("CREATE_INDEX", "True", sep = "=")
        params.set_paths(OUT_BAM = output_bam,
                         OUT_BAI = swap_ext(output_bam, ".bai"))

        for (index, filename) in enumerate(input_bams, start = 1):
            params.push_parameter("I", "%%(IN_BAM_%02i)s" % index, sep = "=")
            params.set_paths("IN_BAM_%02i" % index, filename)
            
        params.set_parameter("SO", "coordinate", sep = "=", fixed = False)

        return {"command" : params}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        description =  "<Merge BAMs: %i file(s) -> '%s'>" \
            % (len(parameters.input_bams), parameters.output_bam)
        CommandNode.__init__(self, 
                             command      = parameters.command.finalize(),
                             description  = description,
                             dependencies = parameters.dependencies)


def concatenate_input_bams(config, input_bams, out = AtomicCmd.PIPE):
    """Transparent concatenation of input BAMs.

    Return a tuple containing a list of nodes (0 or 1), and an
    object which may be passed to the IN_STDIN of an AtomicCmd
    (either an AtomicCmd, or a filename). This allows transparent
    concatenation when multiple files are specified, while
    avoiding needless overhead when there is only 1 input file."""

    input_bams = safe_coerce_to_tuple(input_bams)
    if len(input_bams) == 1:
        return [], input_bams[0]

    jar_file = os.path.join(config.jar_root, "MergeSamFiles.jar")
    params = AtomicJavaParams(config, jar_file)

    if out == AtomicCmd.PIPE:
        params.set_paths(OUT_STDOUT = out)
        params.set_parameter("OUTPUT", "/dev/stdout", sep = "=")
    else:
        params.set_parameter("OUTPUT", out, sep = "=")

    params.set_parameter("CREATE_INDEX", "False", sep = "=")
    params.set_parameter("COMPRESSION_LEVEL",  0, sep = "=")

    for (index, filename) in enumerate(safe_coerce_to_tuple(input_bams), start = 1):
        params.push_parameter("I", "%%(IN_BAM_%02i)s" % index, sep = "=")
        params.set_paths("IN_BAM_%02i" % index, filename)

    params.set_parameter("SO", "coordinate", sep = "=", fixed = False)

    cmd = params.finalize()
    return [cmd], cmd
