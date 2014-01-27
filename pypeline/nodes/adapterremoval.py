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
import itertools
import subprocess

from pypeline.node import \
     CommandNode, \
     NodeError
from pypeline.atomiccmd.command import \
     AtomicCmd,\
     CmdError
from pypeline.atomiccmd.sets import \
     ParallelCmds
from pypeline.atomiccmd.builder import \
     AtomicCmdBuilder, \
     use_customizable_cli_parameters, \
     create_customizable_cli_parameters


import pypeline.common.utilities as utilities
import pypeline.common.fileutils as fileutils
import pypeline.common.versions as versions
import pypeline.common.formats.fastq as fastq
import pypeline.common.sampling as sampling



VERSION_14 = "1.4"
_VERSION_14_CHECK = versions.Requirement(call   = ("AdapterRemoval", "--version"),
                                         search = r"ver. (\d+)\.(\d+)",
                                         checks = versions.EQ(1, 4))

VERSION_15 = "1.5+"
_VERSION_15_CHECK = versions.Requirement(call   = ("AdapterRemoval", "--version"),
                                         search = r"ver. (\d+)\.(\d+)",
                                         checks = versions.EQ(1, 5))


class SE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_files, output_prefix, output_format = "bz2",
                  quality_offset = 33, version = VERSION_15, dependencies = ()):
        # See below for parameters in common between SE/PE
        cmd = _get_common_parameters(version)

        # Uncompressed reads (piped from unicat)
        cmd.set_option("--file1",    "%(TEMP_IN_READS)s")

        # Prefix for output files
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        # Quality offset for Phred (or similar) scores
        assert quality_offset in (33, 64)
        cmd.set_option("--qualitybase", quality_offset)

        basename = os.path.basename(output_prefix)
        cmd.set_kwargs(# Only settings file is saved, rest is temporary files
                       OUT_SETTINGS        = output_prefix + ".settings",
                       TEMP_OUT_BASENAME   = basename,

                       # Named pipe for uncompressed input (unicat)
                       TEMP_IN_READS       = "uncompressed_input",

                       # Named pipes for output of AdapterRemova
                       TEMP_OUT_LINK_1     = basename + ".truncated",
                       TEMP_OUT_LINK_2     = basename + ".discarded")

        return {"basename"      : basename,
                "version"       : version,
                "command"       : cmd}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._quality_offset = parameters.quality_offset
        self._basename = parameters.basename

        zcat           = _build_unicat_command(parameters.input_files, "uncompressed_input")
        zip_truncated  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".truncated")
        zip_discarded  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".discarded")
        adapterrm      = parameters.command.finalize()

        # Opening of pipes block, so the order of these commands is dependent upon
        # the order of file-opens in atomiccmd and the the programs themselves.
        commands = ParallelCmds([adapterrm, zip_discarded, zip_truncated, zcat])
        CommandNode.__init__(self,
                             command      = commands,
                             description  = "<AdapterRM (SE): %s -> '%s.*'>" \
                                 % (fileutils.describe_files(parameters.input_files),
                                    parameters.output_prefix),
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        _check_qualities(self.input_files, self._quality_offset)

        os.mkfifo(os.path.join(temp, self._basename + ".truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, "uncompressed_input"))

        CommandNode._setup(self, config, temp)




class PE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_files_1, input_files_2, output_prefix,
                  quality_offset = 33, output_format = "bz2",
                  version = VERSION_15, dependencies = ()):
        cmd = _get_common_parameters(version)

        # Uncompressed mate 1 and 2 reads (piped from unicat)
        cmd.set_option("--file1",    "%(TEMP_IN_READS_1)s")
        cmd.set_option("--file2",    "%(TEMP_IN_READS_2)s")

        # Prefix for output files, ensure that all end up in temp folder
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        # Quality offset for Phred (or similar) scores
        assert quality_offset in (33, 64)
        cmd.set_option("--qualitybase", quality_offset)

        # Output files are explicity specified, to ensure that the order is the same here
        # as below. A difference in the order in which files are opened can cause a deadlock,
        # due to the use of named pipes (see __init__).
        cmd.set_option("--output1", "%(TEMP_OUT_LINK_PAIR1)s")
        cmd.set_option("--output2", "%(TEMP_OUT_LINK_PAIR2)s")
        cmd.set_option("--outputcollapsed", "%(TEMP_OUT_LINK_ALN)s")
        if version == VERSION_15:
            cmd.set_option("--outputcollapsedtruncated", "%(TEMP_OUT_LINK_ALN_TRUNC)s")
        cmd.set_option("--singleton", "%(TEMP_OUT_LINK_UNALN)s")
        cmd.set_option("--discarded", "%(TEMP_OUT_LINK_DISC)s")

        basename = os.path.basename(output_prefix)
        cmd.set_kwargs(# Only settings file is saved, rest is temporary files
                       OUT_SETTINGS        = output_prefix + ".settings",
                       TEMP_OUT_BASENAME   = basename,

                       # Named pipes for uncompressed input (unicat)
                       TEMP_IN_READS_1     = "uncompressed_input_1",
                       TEMP_IN_READS_2     = "uncompressed_input_2",

                       # Named pipes for output of AdapterRemoval
                       TEMP_OUT_LINK_PAIR1 = basename + ".pair1.truncated",
                       TEMP_OUT_LINK_PAIR2 = basename + ".pair2.truncated",
                       TEMP_OUT_LINK_DISC  = basename + ".discarded")

        if version == VERSION_15:
            cmd.set_kwargs(TEMP_OUT_LINK_ALN       = basename + ".collapsed",
                           TEMP_OUT_LINK_ALN_TRUNC = basename + ".collapsed.truncated",
                           TEMP_OUT_LINK_UNALN     = basename + ".singleton.truncated")
        elif version == VERSION_14:
            cmd.set_kwargs(TEMP_OUT_LINK_ALN       = basename + ".singleton.aln.truncated",
                           TEMP_OUT_LINK_UNALN     = basename + ".singleton.unaln.truncated")
        else:
            assert False

        return {"basename"       : basename,
                "command"        : cmd}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._quality_offset = parameters.quality_offset
        self._version    = parameters.version
        self._basename   = parameters.basename
        if len(parameters.input_files_1) != len(parameters.input_files_2):
            raise CmdError("Number of mate 1 files differ from mate 2 files: %i != %i" \
                               % (len(parameters.input_files_1),
                                  len(parameters.input_files_2)))

        zcat_pair_1    = _build_unicat_command(parameters.input_files_1, "uncompressed_input_1")
        zcat_pair_2    = _build_unicat_command(parameters.input_files_2, "uncompressed_input_2")
        zip_pair_1     = _build_zip_command(parameters.output_format, parameters.output_prefix, ".pair1.truncated")
        zip_pair_2     = _build_zip_command(parameters.output_format, parameters.output_prefix, ".pair2.truncated")
        zip_discarded  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".discarded")
        adapterrm      = parameters.command.finalize()

        commands = [adapterrm, zip_pair_1, zip_pair_2]
        if parameters.version == VERSION_15:
            zip_aln        = _build_zip_command(parameters.output_format, parameters.output_prefix, ".collapsed")
            zip_aln_trunc  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".collapsed.truncated")
            zip_unaligned  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".singleton.truncated")
            commands      += [zip_aln, zip_aln_trunc, zip_unaligned]
        else:
            zip_aln        = _build_zip_command(parameters.output_format, parameters.output_prefix, ".singleton.aln.truncated")
            zip_unaligned  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".singleton.unaln.truncated")
            commands      += [zip_aln, zip_unaligned]
        commands += [zip_discarded, zcat_pair_1, zcat_pair_2]

        # Opening of pipes block, so the order of these commands is dependent upon
        # the order of file-opens in atomiccmd and the the programs themselves.
        commands = ParallelCmds(commands)

        description  = "<AdapterRM (PE): %s -> '%s.*'>" \
            % (fileutils.describe_paired_files(parameters.input_files_1,
                                               parameters.input_files_2),
               parameters.output_prefix)

        CommandNode.__init__(self,
                             command      = commands,
                             description  = description,
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        _check_qualities(self.input_files, self._quality_offset)

        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, self._basename + ".pair1.truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".pair2.truncated"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_1"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_2"))

        if self._version == VERSION_15:
            os.mkfifo(os.path.join(temp, self._basename + ".collapsed"))
            os.mkfifo(os.path.join(temp, self._basename + ".collapsed.truncated"))
            os.mkfifo(os.path.join(temp, self._basename + ".singleton.truncated"))
        else:
            os.mkfifo(os.path.join(temp, self._basename + ".singleton.aln.truncated"))
            os.mkfifo(os.path.join(temp, self._basename + ".singleton.unaln.truncated"))

        CommandNode._setup(self, config, temp)



def _build_unicat_command(input_files, output_file):
    paths = {"TEMP_OUT_CAT" : output_file}
    call = ["unicat", "--output", "%(TEMP_OUT_CAT)s"]
    for (index, filename) in enumerate(utilities.safe_coerce_to_tuple(input_files)):
        key = "IN_CAT_%02i" % index

        call.append("%%(%s)s" % key)
        paths[key] = filename

    return AtomicCmd(call, **paths)


def _build_zip_command(output_format, prefix, name, output = None):
    if output_format == "bz2":
        command, ext = "bzip2", ".bz2"
    elif output_format == "gz":
        command, ext = "gzip", ".gz"
    else:
        raise CmdError("Invalid output-format (%s), please select 'gz' or 'bz2'" \
                       % repr(output_format))

    basename = os.path.basename(prefix)
    return AtomicCmd([command, "-c"],
                     TEMP_IN_STDIN = basename + name,
                     OUT_STDOUT    = prefix + (output or name) + ext)


def _get_common_parameters(version):
    if version == VERSION_14:
        version_check = _VERSION_14_CHECK
    elif version == VERSION_15:
        version_check = _VERSION_15_CHECK
    else:
        raise CmdError("Unknown version: %s" % version)

    cmd = AtomicCmdBuilder("AdapterRemoval",
                           CHECK_VERSION = version_check)

    # Trim Ns at read ends
    cmd.set_option("--trimns", fixed = False)
    # Trim low quality scores
    cmd.set_option("--trimqualities", fixed = False)
    # Offset of quality scores
    cmd.set_option("--qualitybase", 33, fixed = False)
    # Merge pairs where the sequence is overlapping
    cmd.set_option("--collapse", fixed = False)

    return cmd


def _check_qualities(filenames, required_offset):
    for filename in filenames:
        qualities = _read_sequences(filename)
        offsets = fastq.classify_quality_strings(qualities)
        if offsets == fastq.OFFSET_BOTH:
            raise NodeError("FASTQ file contains quality scores with both "
                            "quality offsets (33 and 64); file may be "
                            "unexpected format or corrupt: %r" % (filename,))
        elif offsets == fastq.OFFSET_MISSING:
            raise NodeError("FASTQ file did not contain quality scores; file "
                            "may be unexpected format or corrupt: %r"
                            % (filename,))
        elif offsets not in (fastq.OFFSET_AMBIGIOUS, required_offset):
            raise NodeError("FASTQ file contains quality scores with wrong "
                            "quality score offset (%i); expected reads with "
                            "quality score offset %i: %s"
                            % (offsets, required_offset, filename))


def _read_sequences(filename):
    try:
        unicat = None
        unicat = subprocess.Popen(['unicat', filename],
                                  stderr=subprocess.PIPE,
                                  stdout=subprocess.PIPE)
        qualities = itertools.islice(unicat.stdout, 3, None, 4)

        return sampling.reservoir_sample(qualities, 100000)
    finally:
        # Check return-codes
        rc_unicat = unicat.wait() if unicat else False
        if rc_unicat:
            message = "Error running unicat:\n" \
                      "  Unicat return-code = %i\n\n%s" \
                      % (rc_unicat, unicat.stderr.read())
            raise NodeError(message)
