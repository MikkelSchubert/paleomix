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

import pypeline.tools.unicat as unicat




class SE_AdapterRemovalNode(CommandNode):
    def __init__(self, input_files, output_prefix, dependencies = ()):
        zcat = unicat.build_atomiccmd(AtomicCmd, input_files, "uncompressed_input")
        gzip_truncated = _build_gzip_command(output_prefix, ".truncated")
        gzip_discarded = _build_gzip_command(output_prefix, ".discarded")

        call = ["AdapterRemoval",
                "--mm", 3,
                "--pcr2", "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
                "--minlength", 25,
                "--trimns", 
                "--trimqualities",
                "--qualitybase", 33,
                "--file1",    "%(TEMP_IN_READS)s",
                "--basename", "%(TEMP_OUT_BASENAME)s"]

        self._basename = os.path.basename(output_prefix)
        adapterrm = AtomicCmd(call,
                              OUT_SETTINGS        = output_prefix + ".settings",
                              TEMP_IN_READS       = "uncompressed_input",
                              TEMP_OUT_BASENAME   = self._basename,
                              
                              # Cleanup of symlinks
                              TEMP_OUT_LINK_1     = self._basename + ".truncated",
                              TEMP_OUT_LINK_2     = self._basename + ".discarded",
                              TEMP_OUT_LINK_3     = "uncompressed_input")


        # Opening of pipes block, so the order of these commands is dependent upon
        # the order of file-opens in atomiccmd and the the programs themselves.
        commands = ParallelCmds([adapterrm, gzip_discarded, gzip_truncated, zcat])

        CommandNode.__init__(self,
                             command      = commands,
                             description  = "<SE_AdapterRM: %s>" \
                                 % (self._desc_files(input_files),),
                             dependencies = dependencies)

    def _setup(self, config, temp):
        os.mkfifo(os.path.join(temp, self._basename + ".truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, "uncompressed_input"))

        CommandNode._setup(self, config, temp)



class PE_AdapterRemovalNode(CommandNode):
    def __init__(self, input_files_1, input_files_2, output_prefix, dependencies = ()):
        if len(input_files_1) != len(input_files_2):
            raise CmdError("Number of mate 1 files differ from mate 2 files: %i != %i" \
                               % (len(input_files_1), len(input_files_2)))

        zcat_pair_1    = unicat.build_atomiccmd(AtomicCmd, input_files_1, "uncompressed_input_1")
        zcat_pair_2    = unicat.build_atomiccmd(AtomicCmd, input_files_2, "uncompressed_input_2")
        gzip_pair_1    = _build_gzip_command(output_prefix, ".pair1.truncated")
        gzip_pair_2    = _build_gzip_command(output_prefix, ".pair2.truncated")
        gzip_aligned   = _build_gzip_command(output_prefix, ".singleton.aln.truncated")
        gzip_unaligned = _build_gzip_command(output_prefix, ".singleton.unaln.truncated")
        gzip_discarded = _build_gzip_command(output_prefix, ".discarded")

        call = ["AdapterRemoval",
                "--mm", 3,
                "--pcr2", "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
                "--minlength", 25,
                "--trimns", 
                "--trimqualities",
                "--qualitybase", 33,
                "--collapse",
                "--file1",    "%(TEMP_IN_READS_1)s",
                "--file2",    "%(TEMP_IN_READS_2)s",
                "--basename", "%(TEMP_OUT_BASENAME)s"]

        self._basename = os.path.basename(output_prefix)
        adapterrm = AtomicCmd(call,
                              OUT_SETTINGS        = output_prefix + ".settings",
                              TEMP_IN_READS_1     = "uncompressed_input_1",
                              TEMP_IN_READS_2     = "uncompressed_input_2",
                              TEMP_OUT_BASENAME   = self._basename,
                              
                              # Cleanup of symlinks
                              TEMP_OUT_LINK_1     = self._basename + ".singleton.aln.truncated",
                              TEMP_OUT_LINK_2     = self._basename + ".singleton.unaln.truncated",
                              TEMP_OUT_LINK_3     = self._basename + ".pair1.truncated",
                              TEMP_OUT_LINK_4     = self._basename + ".pair2.truncated",
                              TEMP_OUT_LINK_5     = self._basename + ".discarded",
                              TEMP_OUT_LINK_6     = "uncompressed_input_1",
                              TEMP_OUT_LINK_7     = "uncompressed_input_2")

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

        CommandNode.__init__(self,
                             command      = commands,
                             description  = "<PE_AdapterRM: %s>" \
                                     % (self._desc_files(input_files_1).replace("file", "pair"),),
                             dependencies = dependencies)

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
    return AtomicCmd(["gzip", "-c", "%(TEMP_IN_GZ)s"],
                     TEMP_IN_GZ = basename + name,
                     OUT_STDOUT = prefix + name + ".gz")

