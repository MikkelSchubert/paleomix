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
import shutil
import getpass

from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.builder import \
     AtomicJavaCmdBuilder, \
     create_customizable_cli_parameters, \
     use_customizable_cli_parameters
from pypeline.common.fileutils import \
     swap_ext, \
     try_rmtree, \
     describe_files
from pypeline.common.utilities import safe_coerce_to_tuple
import pypeline.common.versions as versions


_PICARD_VERSION_CACHE = {}
def _picard_version(config, jar_file):
    if jar_file not in _PICARD_VERSION_CACHE:
        params = AtomicJavaCmdBuilder(jar_file,
                                      temp_root = config.temp_root,
                                      jre_options = config.jre_options)
        params.add_value("--version")
        requirement = versions.Requirement(call   = params.finalized_call,
                                           name   = "Picard " + os.path.basename(jar_file),
                                           search = r"^(\d+)\.(\d+)",
                                           checks = versions.GE(1, 82))
        _PICARD_VERSION_CACHE[jar_file] = requirement
    return _PICARD_VERSION_CACHE[jar_file]


class PicardNode(CommandNode):
    """Base class for nodes using Picard Tools; adds an additional cleanup step,
    in order to allow the jars to be run using the same temporary folder as any
    other commands assosiated with the node."""

    def _teardown(self, config, temp):
        # Picard creates a folder named after the user in the temp-root
        try_rmtree(os.path.join(temp, getpass.getuser()))

        CommandNode._teardown(self, config, temp)


class ValidateBAMNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bam, output_log = None, dependencies = ()):
        jar_file = os.path.join(config.jar_root, "ValidateSamFile.jar")
        params = AtomicJavaCmdBuilder(jar_file, jre_options = config.jre_options)

        params.set_option("I", "%(IN_BAM)s", sep = "=")
        params.set_kwargs(IN_BAM     = input_bam,
                         OUT_STDOUT = output_log or swap_ext(input_bam, ".validated"),
                         CHECK_JAR  = _picard_version(config, jar_file))

        return {"command"      : params,
                "dependencies" : dependencies}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        PicardNode.__init__(self,
                            command      = parameters.command.finalize(),
                            description  = "<Validate BAM: '%s'>" % (parameters.input_bam,),
                            dependencies = parameters.dependencies)




class BuildSequenceDictNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, reference, dependencies = ()):
        jar_file = os.path.join(config.jar_root, "CreateSequenceDictionary.jar")
        params = AtomicJavaCmdBuilder(jar_file, jre_options = config.jre_options)

        params.set_option("R", "%(IN_REF)s", sep = "=")
        params.set_option("O", "%(OUT_DICT)s", sep = "=")
        params.set_kwargs(IN_REF     = reference,
                          OUT_DICT   = swap_ext(reference, ".dict"),
                          CHECK_JAR  = _picard_version(config, jar_file))

        return {"command"      : params,
                "dependencies" : dependencies}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        PicardNode.__init__(self,
                            command      = parameters.command.finalize(),
                            description  = "<SequenceDictionary: '%s'>" % (parameters.reference,),
                            dependencies = parameters.dependencies)




class MarkDuplicatesNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bams, output_bam, output_metrics = None, keep_dupes = False, dependencies = ()):
        jar_file = os.path.join(config.jar_root, "MarkDuplicates.jar")
        params = AtomicJavaCmdBuilder(jar_file, jre_options = config.jre_options)

        # Create .bai index, since it is required by a lot of other programs
        params.set_option("CREATE_INDEX", "True", sep = "=")

        params.set_option("OUTPUT", "%(OUT_BAM)s", sep = "=")
        params.set_option("METRICS_FILE", "%(OUT_METRICS)s", sep = "=")

        input_bams = safe_coerce_to_tuple(input_bams)
        for (index, filename) in enumerate(input_bams):
            params.add_option("I", "%%(IN_BAM_%02i)s" % index, sep = "=")
            params.set_kwargs(**{("IN_BAM_%02i" % index) : filename})

        if not keep_dupes:
            # Remove duplicates from output by default to save disk-space
            params.set_option("REMOVE_DUPLICATES", "True", sep = "=", fixed = False)

        params.set_kwargs(OUT_BAM     = output_bam,
                         OUT_BAI     = swap_ext(output_bam, ".bai"),
                         OUT_METRICS = output_metrics or swap_ext(output_bam, ".metrics"),
                         CHECK_JAR  = _picard_version(config, jar_file))

        return {"command"      : params,
                "dependencies" : dependencies}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        description =  "<MarkDuplicates: %s>" % (describe_files(parameters.input_bams),)
        PicardNode.__init__(self,
                            command      = parameters.command.finalize(),
                            description  = description,
                            dependencies = parameters.dependencies)




class MergeSamFilesNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bams, output_bam, dependencies = ()):
        jar_file = os.path.join(config.jar_root, "MergeSamFiles.jar")
        params = AtomicJavaCmdBuilder(jar_file, jre_options = config.jre_options)

        params.set_option("OUTPUT", "%(OUT_BAM)s", sep = "=")
        params.set_option("CREATE_INDEX", "True", sep = "=")
        params.set_kwargs(OUT_BAM = output_bam,
                         OUT_BAI = swap_ext(output_bam, ".bai"),
                         CHECK_JAR  = _picard_version(config, jar_file))

        for (index, filename) in enumerate(input_bams, start = 1):
            params.add_option("I", "%%(IN_BAM_%02i)s" % index, sep = "=")
            params.set_kwargs(**{("IN_BAM_%02i" % index) : filename})

        params.set_option("SO", "coordinate", sep = "=", fixed = False)

        return {"command"      : params,
                "dependencies" : dependencies}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        description =  "<Merge BAMs: %i file(s) -> '%s'>" \
            % (len(parameters.input_bams), parameters.output_bam)
        PicardNode.__init__(self,
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
    params = AtomicJavaCmdBuilder(jar          = jar_file,
                                  temp_root    = config.temp_root,
                                  jre_options  = config.jre_options)
    params.set_kwargs(CHECK_JAR  = _picard_version(config, jar_file))

    if out == AtomicCmd.PIPE:
        params.set_kwargs(OUT_STDOUT = out)
        params.set_option("OUTPUT", "/dev/stdout", sep = "=")
    else:
        params.set_option("OUTPUT", out, sep = "=")

    params.set_option("CREATE_INDEX", "False", sep = "=")
    params.set_option("COMPRESSION_LEVEL",  0, sep = "=")

    for (index, filename) in enumerate(safe_coerce_to_tuple(input_bams), start = 1):
        params.add_option("I", "%%(IN_BAM_%02i)s" % index, sep = "=")
        params.set_kwargs(**{"IN_BAM_%02i" % index : filename})

    params.set_option("SO", "coordinate", sep = "=", fixed = False)

    cmd = params.finalize()
    return [cmd], cmd
