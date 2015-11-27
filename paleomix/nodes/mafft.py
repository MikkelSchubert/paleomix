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
from paleomix.node import \
    CommandNode, \
    NodeError
from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder, \
    use_customizable_cli_parameters, \
    create_customizable_cli_parameters
from paleomix.common.fileutils import \
    reroot_path
from paleomix.common.formats.msa import \
    MSA, \
    MSAError
import paleomix.common.versions as versions


MAFFT_VERSION = versions.Requirement(call   = ("mafft", "--version"),
                                     search = r"v(\d+)\.(\d+)",
                                     checks = versions.GE(7, 0))


# Presets mainly taken from
# http://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html
_PRESETS = {
    "mafft"    : ["mafft"],
    "auto"     : ["mafft", "--auto"],
    "fft-ns-1" : ["mafft-fftns", "--retree", 1],
    "fft-ns-2" : ["mafft-fftns"],
    "fft-ns-i" : ["mafft-fftnsi"],
    "nw-ns-i"  : ["mafft-nwnsi"],
    "l-ins-i"  : ["mafft-linsi"],
    "e-ins-i"  : ["mafft-einsi"],
    "g-ins-i"  : ["mafft-ginsi"],
    }




class MAFFTNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file, output_file, algorithm = "auto", dependencies = ()):
        command = AtomicCmdBuilder(_PRESETS[algorithm.lower()])
        command.add_value("%(IN_FASTA)s")
        command.set_kwargs(IN_FASTA   = input_file,
                           OUT_STDOUT = output_file,
                           CHECK_VERSION = MAFFT_VERSION)

        return {"command"      : command,
                "dependencies" : dependencies}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._output_file = parameters.output_file
        description = "<MAFFTNode (%s): '%s' -> '%s'>" \
                % (parameters.algorithm,
                   parameters.input_file,
                   parameters.output_file)

        CommandNode.__init__(self,
                             command      = parameters.command.finalize(),
                             description  = description,
                             dependencies = parameters.dependencies)

    def _teardown(self, config, temp):
        # Validate output from MAFFT
        output_file = reroot_path(temp, self._output_file)
        try:
            MSA.from_file(output_file)
        except MSAError, error:
            raise NodeError("Invalid MSA produced by MAFFT:\n%s" % (error,))
        CommandNode._teardown(self, config, temp)
