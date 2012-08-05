#!/usr/bin/python

import os

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds





class SE_AdapterRemovalNode(CommandNode):
    def __init__(self, input_files, output_prefix, dependencies = ()):
        zcat = _build_zcat_command(input_files, "uncompressed_input")
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
                             description  = "<SE_AdapterRM: %s -> '%s.*'>" \
                                     % (_desc_input(input_files), output_prefix),
                             dependencies = dependencies)

    def _setup(self, config, temp):
        os.mkfifo(os.path.join(temp, self._basename + ".truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, "uncompressed_input"))

        CommandNode._setup(self, config, temp)



class PE_AdapterRemovalNode(CommandNode):
    def __init__(self, input_files_1, input_files_2, output_prefix, dependencies = ()):
        zcat_pair_1    = _build_zcat_command(input_files_1, "uncompressed_input_1")
        zcat_pair_2    = _build_zcat_command(input_files_2, "uncompressed_input_2")
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
                             description  = "<PE_AdapterRM: %s -> '%s.*'>" \
                                     % (_desc_input(input_files_1), output_prefix),
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
#    return AtomicCmd(["gzip", "-c", "%(TEMP_IN)s"],
#                     TEMP_IN    = basename + name,
#                     OUT_STDOUT = prefix + name + ".gz")
    return AtomicCmd(["gzip"],
                     TEMP_IN_STDIN = basename + name,
                     OUT_STDOUT    = prefix + name + ".gz")


def _build_zcat_command(input_files, output_file):
        counter = 1
        zcat_call = ["zcat"]
        zcat_dict = {"TEMP_OUT_STDOUT" : output_file}
        for filename in sorted(input_files):
            zcat_call.append(filename)
            zcat_dict["IN_FILE_%04i" % counter] = filename
            counter += 1

        return AtomicCmd(zcat_call, **zcat_dict)


def _desc_input(files):
    paths = set(os.path.dirname(filename) for filename in files)
    if len(paths) == 1:
        return "%i file(s) in '%s'" % (len(files), paths.pop())
    else:
        return "%i file(s)" % (len(files),)
