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
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
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
import getpass

from paleomix.node import CommandNode
from paleomix.atomiccmd.builder import \
    AtomicJavaCmdBuilder, \
    create_customizable_cli_parameters, \
    use_customizable_cli_parameters
from paleomix.common.fileutils import \
    swap_ext, \
    try_rmtree, \
    try_remove, \
    reroot_path, \
    describe_files
from paleomix.common.utilities import \
    safe_coerce_to_tuple
import paleomix.common.versions as versions
import paleomix.common.system


class PicardNode(CommandNode):
    """Base class for nodes using Picard Tools; adds an additional cleanup
    step, in order to allow the jars to be run using the same temporary folder
    as any other commands associated with the node. This is nessesary as some
    Picard tools create a large number of temporary files, leading to potential
    performance issues if these are located in the same folder.
    """

    def _teardown(self, config, temp):
        # Picard creates a folder named after the user in the temp-root
        try_rmtree(os.path.join(temp, getpass.getuser()))
        # Some JREs may create a folder for temporary performance counters
        try_rmtree(os.path.join(temp, "hsperfdata_" + getpass.getuser()))

        CommandNode._teardown(self, config, temp)


class ValidateBAMNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bam, output_log=None, dependencies=()):
        params = picard_command(config, "ValidateSamFile")
        _set_max_open_files(params, "MAX_OPEN_TEMP_FILES")

        params.set_option("I", "%(IN_BAM)s", sep="=")

        output_log = output_log or swap_ext(input_bam, ".validated")
        params.set_kwargs(IN_BAM=input_bam,
                          OUT_STDOUT=output_log)

        return {"command": params,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        description = "<Validate BAM: '%s'>" % (parameters.input_bam,)
        PicardNode.__init__(self,
                            command=parameters.command.finalize(),
                            description=description,
                            dependencies=parameters.dependencies)


class BuildSequenceDictNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, reference, dependencies=()):
        params = picard_command(config, "CreateSequenceDictionary")

        params.set_option("R", "%(TEMP_OUT_REF)s", sep="=")
        params.set_option("O", "%(OUT_DICT)s", sep="=")
        params.set_kwargs(IN_REF=reference,
                          TEMP_OUT_REF=os.path.basename(reference),
                          OUT_DICT=swap_ext(reference, ".dict"))

        return {"command": params,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._in_reference = os.path.abspath(parameters.reference)
        description = "<SequenceDictionary: '%s'>" % (parameters.reference,)

        PicardNode.__init__(self,
                            command=parameters.command.finalize(),
                            description=description,
                            dependencies=parameters.dependencies)

    def _setup(self, _config, temp):
        os.symlink(self._in_reference, reroot_path(temp, self._in_reference))


class MarkDuplicatesNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bams, output_bam, output_metrics=None,
                  keep_dupes=False, dependencies=()):
        params = picard_command(config, "MarkDuplicates")
        _set_max_open_files(params, "MAX_FILE_HANDLES")

        # Create .bai index, since it is required by a lot of other programs
        params.set_option("CREATE_INDEX", "True", sep="=")

        params.set_option("OUTPUT", "%(OUT_BAM)s", sep="=")
        params.set_option("METRICS_FILE", "%(OUT_METRICS)s", sep="=")
        params.add_multiple_options("I", input_bams, sep="=")

        if not keep_dupes:
            # Remove duplicates from output by default to save disk-space
            params.set_option("REMOVE_DUPLICATES", "True",
                              sep="=", fixed=False)

        output_metrics = output_metrics or swap_ext(output_bam, ".metrics")
        params.set_kwargs(OUT_BAM=output_bam,
                          OUT_BAI=swap_ext(output_bam, ".bai"),
                          OUT_METRICS=output_metrics)

        return {"command": params,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        description = "<MarkDuplicates: %s>" \
            % (describe_files(parameters.input_bams),)
        PicardNode.__init__(self,
                            command=parameters.command.finalize(),
                            description=description,
                            dependencies=parameters.dependencies)


class MergeSamFilesNode(PicardNode):
    @create_customizable_cli_parameters
    def customize(cls, config, input_bams, output_bam, dependencies=()):
        params = picard_command(config, "MergeSamFiles")

        params.set_option("OUTPUT", "%(OUT_BAM)s", sep="=")
        params.set_option("CREATE_INDEX", "True", sep="=")
        params.set_option("SO", "coordinate", sep="=", fixed=False)
        params.add_multiple_options("I", input_bams, sep="=")

        params.set_kwargs(OUT_BAM=output_bam,
                          OUT_BAI=swap_ext(output_bam, ".bai"))

        return {"command": params,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        description = "<Merge BAMs: %i file(s) -> '%s'>" \
            % (len(parameters.input_bams), parameters.output_bam)
        PicardNode.__init__(self,
                            command=parameters.command.finalize(),
                            description=description,
                            dependencies=parameters.dependencies)


class MultiBAMInput(object):
    """Container used to ease processing of 1 or more BAM files; used in
    conjunctin with MultiBAMInputNode.
    """

    def __init__(self, config, input_bams, pipename="input.bam", indexed=True):
        self.pipe = pipename
        self.indexed = indexed
        self.files = safe_coerce_to_tuple(input_bams)

        self.commands = []
        self.kwargs = {"TEMP_IN_BAM": self.pipe}
        if len(self.files) > 1:
            params = picard_command(config, "MergeSamFiles")

            params.set_option("SO", "coordinate", sep="=", fixed=False)
            params.set_option("CREATE_INDEX", "False", sep="=")
            params.set_option("COMPRESSION_LEVEL", 0, sep="=")
            params.set_option("OUTPUT", "%(TEMP_OUT_BAM)s", sep="=")
            params.add_multiple_options("I", input_bams, sep="=")

            params.set_kwargs(TEMP_OUT_BAM=self.pipe)

            self.commands = [params.finalize()]
        else:
            # Ensure that the actual command depends on the input
            self.kwargs["IN_FILE_00"] = self.files[0]

            if indexed:
                self.kwargs["IN_FILE_01"] = swap_ext(self.files[0], ".bai")

    def setup(self, command):
        command.set_kwargs(**self.kwargs)


class MultiBAMInputNode(CommandNode):
    """Node which provides concatenation of input BAM files. Takes a
    MultiBAMInput object, and creates a pipe in the temporary folder which
    yields the concatenated BAM resulting from the concatenation of all input
    files. To avoid unnessary overhead, a symbolic link is used in the case
    where there is only a single input file.

    Usage example:
      class ExampleNode(MultiBAMInputNode):
        def __init__(self, config, input_bams):
            bam_input = MultiBAMInput(config, input_bams)
            command = AtomicCmd(['analyse_bam', '%(TEMP_IN_BAM)s'],
                                TEMP_IN_BAM=bam_input.pipe)
            commands = ParallelCmds(bam_input.commands + [command])
            MultiBAMInputNode.__init__(bam_input=bam_input,
                                       command=commands)
    """

    def __init__(self, bam_input, *args, **kwargs):
        self._bam_input = bam_input
        CommandNode.__init__(self, *args, **kwargs)

    def _setup(self, config, temp_root):
        CommandNode._setup(self, config, temp_root)
        dst_fname = os.path.join(temp_root, self._bam_input.pipe)
        if len(self._bam_input.files) > 1:
            os.mkfifo(dst_fname)
        else:
            src_fname, = self._bam_input.files
            os.symlink(os.path.join(os.getcwd(), src_fname), dst_fname)

            if self._bam_input.indexed:
                src_fname = os.path.join(os.getcwd(), swap_ext(src_fname, ".bai"))
                os.symlink(src_fname, dst_fname + ".bai")

    def _teardown(self, config, temp_root):
        pipe_fname = os.path.join(temp_root, self._bam_input.pipe)
        os.remove(pipe_fname)
        try_remove(pipe_fname + ".bai")
        CommandNode._teardown(self, config, temp_root)


###############################################################################

_PICARD_JAR = "picard.jar"
_PICARD_VERSION_CACHE = {}


def picard_command(config, command):
    """Returns basic AtomicJavaCmdBuilder for Picard tools commands."""
    jar_path = os.path.join(config.jar_root, _PICARD_JAR)

    if jar_path not in _PICARD_VERSION_CACHE:
        params = AtomicJavaCmdBuilder(jar_path,
                                      temp_root=config.temp_root,
                                      jre_options=config.jre_options)

        # Arbitrary command, since just '--version' does not work
        params.set_option("MarkDuplicates")
        params.set_option("--version")

        requirement = versions.Requirement(call=params.finalized_call,
                                           name="Picard tools",
                                           search=r"^(\d+)\.(\d+)",
                                           checks=versions.GE(1, 124))
        _PICARD_VERSION_CACHE[jar_path] = requirement

    version = _PICARD_VERSION_CACHE[jar_path]
    params = AtomicJavaCmdBuilder(jar_path,
                                  temp_root=config.temp_root,
                                  jre_options=config.jre_options,
                                  CHECK_JAR=version)
    params.set_option(command)

    return params


# Fraction of per-process max open files to use
_FRAC_MAX_OPEN_FILES = 0.95


def _set_max_open_files(params, key):
    """Sets the maximum number of open files a picard process
    should use, at most. Conservatively lowered than the actual
    ulimit.
    """
    max_open_files = paleomix.common.system.get_max_open_files()
    if max_open_files is not None:
        max_open_files = int(max_open_files * _FRAC_MAX_OPEN_FILES)
        params.set_option(key, max_open_files, sep="=")
