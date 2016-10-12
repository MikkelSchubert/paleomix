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

from paleomix.node import CommandNode, NodeError
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder, \
    use_customizable_cli_parameters, \
    create_customizable_cli_parameters

from paleomix.atomiccmd.sets import ParallelCmds
from paleomix.nodes.samtools import SAMTOOLS_VERSION
from paleomix.common.fileutils import \
    describe_paired_files, \
    missing_files

import paleomix.common.versions as versions
import paleomix.tools.factory as factory


BWA_VERSION = versions.Requirement(call=("bwa",),
                                   search=r"Version: (\d+)\.(\d+)\.(\d+)",
                                   checks=versions.Or(versions.EQ(0, 5, 9),
                                                      versions.EQ(0, 5, 10),
                                                      versions.EQ(0, 6, 2),
                                                      versions.GE(0, 7, 9)))

BWA_VERSION_07x = versions.Requirement(call=("bwa",),
                                       search=r"Version: (\d+)\.(\d+)\.(\d+)",
                                       checks=versions.GE(0, 7, 9))


class BWAIndexNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file, prefix=None, dependencies=()):
        prefix = prefix if prefix else input_file
        params = _get_bwa_template(("bwa", "index"), prefix, iotype="OUT",
                                   IN_FILE=input_file,
                                   TEMP_OUT_PREFIX=os.path.basename(prefix),
                                   CHECK_BWA=BWA_VERSION)

        # Input fasta sequence
        params.add_value("%(IN_FILE)s")
        # Destination prefix, in temp folder
        params.set_option("-p", "%(TEMP_OUT_PREFIX)s")

        return {"prefix": prefix,
                "command": params,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()
        description = "<BWA Index '%s' -> '%s.*'>" % (parameters.input_file,
                                                      parameters.prefix)
        CommandNode.__init__(self,
                             command=command,
                             description=description,
                             dependencies=parameters.dependencies)


class BWABacktrack(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file, output_file, reference, prefix, threads=1,
                  dependencies=()):
        _check_bwa_prefix(reference)
        threads = _get_max_threads(reference, threads)

        aln_in = _build_cat_command(input_file, "uncompressed_input")
        aln = _get_bwa_template(("bwa", "aln"), prefix,
                                TEMP_IN_FILE="uncompressed_input",
                                OUT_STDOUT=output_file,
                                CHECK_BWA=BWA_VERSION)
        aln.add_value(prefix)
        aln.add_value("%(TEMP_IN_FILE)s")
        aln.set_option("-t", threads)

        return {"commands": {"aln_in": aln_in, "aln": aln},
                "order": ["aln_in", "aln"],
                "threads": threads,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = ParallelCmds([parameters.commands[key].finalize()
                               for key in parameters.order])

        description \
            = _get_node_description(name="BWA",
                                    algorithm='Backtrack',
                                    input_files_1=parameters.input_file,
                                    prefix=parameters.prefix,
                                    threads=parameters.threads)

        CommandNode.__init__(self,
                             command=command,
                             description=description,
                             threads=parameters.threads,
                             dependencies=parameters.dependencies)

    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "uncompressed_input"))


class BWASamse(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file_fq, input_file_sai, output_file,
                  reference, prefix, dependencies=()):
        _check_bwa_prefix(reference)

        samse_in = _build_cat_command(input_file_fq, "uncompressed_input")
        samse = _get_bwa_template(("bwa", "samse"), prefix,
                                  IN_FILE_SAI=input_file_sai,
                                  TEMP_IN_FQ="uncompressed_input",
                                  OUT_STDOUT=AtomicCmd.PIPE,
                                  CHECK_BWA=BWA_VERSION)
        samse.add_value(prefix)
        samse.add_value("%(IN_FILE_SAI)s")
        samse.add_value("%(TEMP_IN_FQ)s")

        order, commands = _process_output(samse, output_file, reference)
        commands["sam_in"] = samse_in
        commands["sam"] = samse

        return {"commands": commands,
                "order": ["sam_in", "sam"] + order,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = ParallelCmds([parameters.commands[key].finalize()
                               for key in parameters.order])

        input_file = parameters.input_file_fq
        description = _get_node_description(name="BWA Samse",
                                            input_files_1=input_file,
                                            prefix=parameters.prefix)

        CommandNode.__init__(self,
                             command=command,
                             description=description,
                             dependencies=parameters.dependencies)

    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "uncompressed_input"))


class BWASampe(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls,
                  input_file_fq_1, input_file_fq_2,
                  input_file_sai_1, input_file_sai_2,
                  output_file, reference, prefix, dependencies=()):
        _check_bwa_prefix(reference)

        sampe_in_1 = _build_cat_command(input_file_fq_1,
                                        "uncompressed_input_1")
        sampe_in_2 = _build_cat_command(input_file_fq_2,
                                        "uncompressed_input_2")

        sampe = _get_bwa_template(("bwa", "sampe"), prefix,
                                  IN_FILE_SAI_1=input_file_sai_1,
                                  IN_FILE_SAI_2=input_file_sai_2,
                                  TEMP_IN_FQ_1="uncompressed_input_1",
                                  TEMP_IN_FQ_2="uncompressed_input_2",
                                  OUT_STDOUT=AtomicCmd.PIPE,
                                  CHECK_BWA=BWA_VERSION)
        sampe.add_value(prefix)
        sampe.add_value("%(IN_FILE_SAI_1)s")
        sampe.add_value("%(IN_FILE_SAI_2)s")
        sampe.add_value("%(TEMP_IN_FQ_1)s")
        sampe.add_value("%(TEMP_IN_FQ_2)s")

        order, commands = _process_output(sampe, output_file, reference,
                                          run_fixmate=True)
        commands["sam_in_1"] = sampe_in_1
        commands["sam_in_2"] = sampe_in_2
        commands["sam"] = sampe

        return {"commands": commands,
                "order": ["sam_in_1", "sam_in_2", "sam"] + order,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = ParallelCmds([parameters.commands[key].finalize()
                               for key in parameters.order])

        input_file_1 = parameters.input_file_fq_1
        input_file_2 = parameters.input_file_fq_2
        description = _get_node_description(name="BWA Sampe",
                                            input_files_1=input_file_1,
                                            input_files_2=input_file_2,
                                            prefix=parameters.prefix)

        CommandNode.__init__(self,
                             command=command,
                             description=description,
                             dependencies=parameters.dependencies)

    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "uncompressed_input_1"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_2"))


class BWAAlgorithmNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file_1, output_file, reference, prefix,
                  input_file_2=None, threads=1, algorithm="mem",
                  dependencies=()):
        if algorithm not in ("mem", "bwasw"):
            raise NotImplementedError("BWA algorithm %r not implemented"
                                      % (algorithm,))

        threads = _get_max_threads(reference, threads)

        zcat_1 = _build_cat_command(input_file_1, "uncompressed_input_1")
        aln = _get_bwa_template(("bwa", algorithm), prefix,
                                TEMP_IN_FILE_1="uncompressed_input_1",
                                OUT_STDOUT=AtomicCmd.PIPE,
                                CHECK_BWA=BWA_VERSION_07x)
        aln.add_value(prefix)
        aln.add_value("%(TEMP_IN_FILE_1)s")

        _, commands = _process_output(aln, output_file, reference)
        commands["aln"] = aln
        commands["zcat_1"] = zcat_1
        if input_file_2:
            aln.add_value("%(TEMP_IN_FILE_2)s")
            aln.set_kwargs(**{"TEMP_IN_FILE_2": "uncompressed_input_2"})
            zcat_2 = _build_cat_command(input_file_2, "uncompressed_input_2")
            commands["zcat_2"] = zcat_2
        else:
            # Ensure that the pipe is automatically removed
            aln.set_kwargs(**{"TEMP_OUT_FILE_2": "uncompressed_input_2"})

        aln.set_option("-t", threads)
        # Mark alternative hits as secondary; required by e.g. Picard
        aln.set_option("-M")

        commands["aln"] = aln
        return {"commands": commands,
                "threads": threads,
                "dependencies": dependencies}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        _check_bwa_prefix(parameters.prefix)
        algorithm = parameters.algorithm.upper()
        algorithm += "_PE" if parameters.input_file_2 else "_SE"
        desc = _get_node_description(name="BWA",
                                     algorithm=algorithm,
                                     input_files_1=parameters.input_file_1,
                                     input_files_2=parameters.input_file_2,
                                     prefix=parameters.prefix)

        command = ParallelCmds([cmd.finalize()
                                for cmd in parameters.commands.itervalues()])
        CommandNode.__init__(self,
                             command=command,
                             description=desc,
                             threads=parameters.threads,
                             dependencies=parameters.dependencies)

    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "uncompressed_input_1"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_2"))


def _process_output(stdin, output_file, reference, run_fixmate=False):
    convert = factory.new("cleanup")
    if reference is not None:
        convert.set_option("--fasta", "%(IN_FASTA_REF)s")
    convert.set_option("--temp-prefix", "%(TEMP_OUT_PREFIX)s")
    convert.set_kwargs(IN_STDIN=stdin,
                       IN_FASTA_REF=reference,
                       OUT_STDOUT=output_file,
                       TEMP_OUT_PREFIX="bam_cleanup",
                       CHECK_SAMTOOLS=SAMTOOLS_VERSION)

    if run_fixmate:
        convert.set_option('--paired-end')

    try:
        if SAMTOOLS_VERSION.version >= (1,):
            convert.set_option('--samtools1x', 'yes')
        else:
            convert.set_option('--samtools1x', 'no')
    except versions.VersionRequirementError:
        pass

    return ["convert"], {"convert": convert}


def _get_bwa_template(call, prefix, iotype="IN", **kwargs):
    extensions = ["amb", "ann", "bwt", "pac", "sa"]
    try:
        if BWA_VERSION.version < (0, 6, 0):
            extensions.extend(("rbwt", "rpac", "rsa"))
    except versions.VersionRequirementError:
        pass  # Ignored here, handled elsewhere

    params = AtomicCmdBuilder(call, **kwargs)
    for postfix in extensions:
        key = "%s_PREFIX_%s" % (iotype, postfix.upper())
        params.set_kwargs(**{key: (prefix + "." + postfix)})

    return params


def _get_max_threads(reference, threads):
    """Returns the maximum number of threads to use when mapping against a
    given reference sequence. This is done since very little gain is obtained
    when using multiple threads for a small genome (e.g. < 1MB). If the
    reference falls below this size, only 1 thread is used (returned),
    otherwise the requested number of threads is returned.
    """
    if reference not in _PREFIX_SIZE_CACHE:
        if reference is None or not os.path.exists(reference):
            _PREFIX_SIZE_CACHE[reference] = None
        else:
            _PREFIX_SIZE_CACHE[reference] = os.path.getsize(reference)

    prefix_size = _PREFIX_SIZE_CACHE[reference]
    if prefix_size is None or prefix_size >= 2 ** 20:  # > 1 MB
        return threads
    return 1
_PREFIX_SIZE_CACHE = {}


def _check_bwa_prefix(prefix):
    """Checks that a given prefix is compatible with the currently
    installed version of BWA. This is required in order to allow
    auto-indexing of prefixes, as indexes produced by v0.5.x and
    by 0.6+ are not only incompatible, but differs in the files
    produced, with 0.5.x producing a handful of additional files.

    As a consequence, simply using normal input-file dependencies
    would result in prefixes being re-indexed if the version of
    BWA was changed from 0.6+ to 0.5.x, and in failures during
    runtime if the version was changed from 0.5.x to 0.6+.

    This function treats that a difference in the version of BWA
    installed and the version implied by the prefix files is an
    error, and therefore requires user intervention."""
    if prefix in _PREFIXES_CHECKED:
        return
    _PREFIXES_CHECKED.add(prefix)

    try:
        bwa_version = BWA_VERSION.version
    except versions.VersionRequirementError:
        return  # Ignored here, reported elsewhere

    # Files unique to v0.5.x
    v05x_files = set((prefix + ext) for ext in (".rbwt", ".rpac", ".rsa"))
    # Files common to v0.5.x, v0.6.x, and v0.7.x
    common_files = set((prefix + ext)
                       for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"))
    all_files = v05x_files | common_files
    current_files = all_files - set(missing_files(all_files))

    expected_version = None
    if (current_files & common_files):
        if bwa_version >= (0, 6, 0):
            if (current_files & v05x_files):
                expected_version = "v0.5.x"
        elif bwa_version < (0, 6, 0):
            if not (current_files & v05x_files):
                expected_version = "v0.6.x or later"

    if expected_version:
        raise NodeError("BWA version is v%s, but prefix appears to be created using %s!\n"
                        "  Your copy of BWA may have changed, or you may be using the wrong\n"
                        "  prefix. To resolve this issue, either change your prefix, re-install\n"
                        "  BWA %s, or remove the prefix files at\n"
                        "    $ ls %s.*" \
                        % (".".join(map(str, bwa_version)), expected_version, expected_version, prefix))
_PREFIXES_CHECKED = set()


def _build_cat_command(input_file, output_file):
    cat = factory.new("cat")
    cat.set_option("--output", "%(TEMP_OUT_CAT)s")
    cat.add_value("%(IN_ARCHIVE)s")
    cat.set_kwargs(TEMP_OUT_CAT=output_file,
                   IN_ARCHIVE=input_file)
    return cat


def _get_node_description(name, input_files_1, input_files_2=None,
                          algorithm=None, prefix=None, threads=1):
    info = []
    if prefix is not None:
        prefix = os.path.basename(prefix)
        if prefix.endswith(".fasta") or prefix.endswith(".fa"):
            prefix = prefix.rsplit(".", 1)[0]

        info.append(prefix)

    if algorithm is not None:
        info.append(algorithm)

    if threads > 1:
        info.append("%i threads" % (threads,))

    file_desc = describe_paired_files(input_files_1, input_files_2 or ())

    return "<%s (%s): %s>" % (name, ", ".join(info), file_desc)
