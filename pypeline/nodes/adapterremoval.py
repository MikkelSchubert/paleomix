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
from pypeline.atomiccmd import AtomicCmd, CmdError
from pypeline.atomicset import ParallelCmds
from pypeline.atomicparams import *
from pypeline.commands import unicat 




class SE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_files, output_prefix, dependencies = ()):
        # See below for parameters in common between SE/PE
        cmd = _get_common_parameters()

        # Uncompressed reads (piped from unicat)
        cmd.set_parameter("--file1",    "%(TEMP_IN_READS)s")

        # Prefix for output files
        cmd.set_parameter("--basename", "%(TEMP_OUT_BASENAME)s")

        basename = os.path.basename(output_prefix)
        cmd.set_paths(# Only settings file is saved, rest is temporary files
                      OUT_SETTINGS        = output_prefix + ".settings",
                      TEMP_OUT_BASENAME   = basename,

                      # Named pipe for uncompressed input (unicat)
                      TEMP_IN_READS       = "uncompressed_input",
                              
                      # Named pipes for output of AdapterRemova
                      TEMP_OUT_LINK_1     = basename + ".truncated",
                      TEMP_OUT_LINK_2     = basename + ".discarded",
                      TEMP_OUT_LINK_3     = "uncompressed_input")

        return {"basename"      : basename,
                "command"       : cmd}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._basename = parameters.basename

        zcat           = unicat.build_atomiccmd(AtomicCmd, parameters.input_files, "uncompressed_input")
        gzip_truncated = _build_gzip_command(parameters.output_prefix, ".truncated")
        gzip_discarded = _build_gzip_command(parameters.output_prefix, ".discarded")
        adapterrm      = parameters.command.finalize()

        # Opening of pipes block, so the order of these commands is dependent upon
        # the order of file-opens in atomiccmd and the the programs themselves.
        commands = ParallelCmds([adapterrm, gzip_discarded, gzip_truncated, zcat])
        CommandNode.__init__(self,
                             command      = commands,
                             description  = "<SE_AdapterRM: %s -> '%s.*'>" \
                                 % (self._desc_files(parameters.input_files),
                                    parameters.output_prefix),
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        os.mkfifo(os.path.join(temp, self._basename + ".truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, "uncompressed_input"))

        CommandNode._setup(self, config, temp)




class PE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(self, input_files_1, input_files_2, output_prefix, dependencies = ()):
        cmd = _get_common_parameters()
        # Merge pairs where the sequence is overlapping
        cmd.set_parameter("--collapse")

        # Uncompressed mate 1 and 2 reads (piped from unicat)
        cmd.set_parameter("--file1",    "%(TEMP_IN_READS_1)s")
        cmd.set_parameter("--file2",    "%(TEMP_IN_READS_2)s")

        # Prefix for output files
        cmd.set_parameter("--basename", "%(TEMP_OUT_BASENAME)s")

        basename = os.path.basename(output_prefix)
        cmd.set_paths(# Only settings file is saved, rest is temporary files
                      OUT_SETTINGS        = output_prefix + ".settings",
                      TEMP_OUT_BASENAME   = basename,
                      
                      # Named pipes for uncompressed input (unicat)
                      TEMP_IN_READS_1     = "uncompressed_input_1",
                      TEMP_IN_READS_2     = "uncompressed_input_2",
                              
                      # Named pipes for output of AdapterRemoval
                      TEMP_OUT_LINK_1     = basename + ".singleton.aln.truncated",
                      TEMP_OUT_LINK_2     = basename + ".singleton.unaln.truncated",
                      TEMP_OUT_LINK_3     = basename + ".pair1.truncated",
                      TEMP_OUT_LINK_4     = basename + ".pair2.truncated",
                      TEMP_OUT_LINK_5     = basename + ".discarded",
                      TEMP_OUT_LINK_6     = "uncompressed_input_1",
                      TEMP_OUT_LINK_7     = "uncompressed_input_2")

        return {"basename"       : basename,
                "command"        : cmd}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._basename = parameters.basename
        if len(parameters.input_files_1) != len(parameters.input_files_2):
            raise CmdError("Number of mate 1 files differ from mate 2 files: %i != %i" \
                               % (len(parameters.input_files_1), 
                                  len(parameters.input_files_2)))

        zcat_pair_1    = unicat.build_atomiccmd(AtomicCmd, parameters.input_files_1, "uncompressed_input_1")
        zcat_pair_2    = unicat.build_atomiccmd(AtomicCmd, parameters.input_files_2, "uncompressed_input_2")
        gzip_pair_1    = _build_gzip_command(parameters.output_prefix, ".pair1.truncated")
        gzip_pair_2    = _build_gzip_command(parameters.output_prefix, ".pair2.truncated")
        gzip_aligned   = _build_gzip_command(parameters.output_prefix, ".singleton.aln.truncated")
        gzip_unaligned = _build_gzip_command(parameters.output_prefix, ".singleton.unaln.truncated")
        gzip_discarded = _build_gzip_command(parameters.output_prefix, ".discarded")
        adapterrm      = parameters.command.finalize()
        
        # Opening of pipes block, so the order of these commands is dependent upon
        # the order of file-opens in atomiccmd and the the programs themselves. 
        commands = ParallelCmds([adapterrm, 
                                 gzip_discarded,
                                 gzip_pair_1,
                                 gzip_pair_2,
                                 gzip_aligned,
                                 gzip_unaligned,
                                 zcat_pair_1,
                                 zcat_pair_2])

        description  = "<PE_AdapterRM: %s -> '%s.*'>" \
            % (self._desc_files(parameters.input_files_1).replace("file", "pair"),
               parameters.output_prefix)

        CommandNode.__init__(self,
                             command      = commands,
                             description  = description,
                             dependencies = parameters.dependencies)


    def _setup(self, config, temp):
        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, self._basename + ".pair1.truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".pair2.truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".singleton.aln.truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".singleton.unaln.truncated"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_1"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_2"))

        CommandNode._setup(self, config, temp)




def _build_gzip_command(prefix, name):
    basename = os.path.basename(prefix)
    return AtomicCmd(["gzip", "-c", "-n"],
                     TEMP_IN_STDIN = basename + name,
                     OUT_STDOUT    = prefix + name + ".gz")


def _get_common_parameters():
    cmd = AtomicParams("AdapterRemoval")

    # Allow 1/3 mismatches in the aligned region
    cmd.set_parameter("--mm", 3, fixed = False)
    # Reverse complement of adapter (required for SE?) 
    cmd.set_parameter("--pcr2", "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", fixed = False)
    # Minimum length of trimmed reads
    cmd.set_parameter("--minlength", 25, fixed = False)
    # Trim Ns at read ends
    cmd.set_parameter("--trimns", fixed = False)
    # Trim low quality scores
    cmd.set_parameter("--trimqualities", fixed = False)
    # Offset of quality scores
    cmd.set_parameter("--qualitybase", 33, fixed = False)

    return cmd
