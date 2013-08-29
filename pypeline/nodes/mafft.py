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

from pypeline.node import CommandNode, MetaNode
from pypeline.atomiccmd.command import AtomicCmd


# Presets mainly taken from
# http://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html
_PRESETS = {
    "auto"     : ["mafft", "--auto"],
    "FFT-NS-1" : ["fftns", "--retree", 1],
    "FFT-NS-2" : ["fftns"],
    "FFT-NS-i" : ["fftnsi", "--maxiterate", 1000],
    "NW-NS-i"  : ["nwnsi",  "--maxiterate", 1000],
    "L-INS-i"  : ["linsi",  "--maxiterate", 1000],
    "E-INS-i"  : ["einsi",  "--maxiterate", 1000],
    "G-INS-i"  : ["ginsi",  "--maxiterate", 1000]}



class MAFFTNode(CommandNode):
    def __init__(self, infile, outfile, preset = "auto", dependencies = ()):
        call = _PRESETS[preset] + ["--quiet", "%(IN_FASTA)s"]
        command = AtomicCmd(call,
                            IN_FASTA   = infile,
                            OUT_STDOUT = outfile)

        description = "<MAFFTNode (%s): '%s' -> '%s'>" \
                % (preset, infile, outfile)

        CommandNode.__init__(self,
                             command      = command,
                             description  = description,
                             dependencies = dependencies)


class MetaMAFFTNode(MetaNode):
    def __init__(self, rootdir, sequences, preset = "auto", subnodes = (), dependencies = ()):
        self._nodemap = {}
        for sequence in sequences:
            prefix      = os.path.join(rootdir, sequence)
            input_file  = prefix + ".fasta"
            output_file = prefix + ".afa"

            self._nodemap[output_file] \
              = MAFFTNode(infile       = input_file,
                          outfile      = output_file,
                          preset       = preset,
                          dependencies = dependencies)

        MetaNode.__init__(self,
                          description  = "<MAFFTAlignSequences (%s): In '%s'>" \
                              % (preset, rootdir),
                          subnodes     = self._nodemap.values(),
                          dependencies = dependencies)

    @property
    def files_to_nodes_map(self):
        return dict(self._nodemap)
