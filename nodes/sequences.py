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
from __future__ import with_statement

import os
import textwrap

import pysam
import pypeline.common.fileutils as fileutils

from pypeline.node import Node



class CollectSequencesNode(Node):
    def __init__(self, infiles, destination, sequences, dependencies):
        self._infiles     = dict(infiles)
        self._destination = str(destination)
        self._sequences   = list(sequences)
        self._outfiles    = []

        for sequence in self._sequences:
            self._outfiles.append(os.path.join(destination, sequence + ".fasta"))

        Node.__init__(self,
                      description  = "<CollectSequences: %i sequences from %i files -> '%s'>" \
                            % (len(sequences), len(self._infiles), destination),
                      input_files  = self._infiles.values(),
                      output_files = self._outfiles,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        handles = []
        for (name, source) in sorted(self._infiles.items()):
            handles.append((name, pysam.Fastafile(source)))

        for sequence in self._sequences:
            filename = os.path.join(temp, sequence + ".fasta")

            with open(filename, "w") as fasta:
                for (name, handle) in handles:
                    fastaseq = textwrap.fill(handle.fetch(sequence), 60)
                    assert fastaseq, (name, sequence)
                    fasta.write(">%s\n%s\n" % (name, fastaseq))

        for (name, handle) in handles:
            handle.close()


    def _teardown(self, _config, temp):
        for sequence in self._sequences:
            filename = sequence + ".fasta"
            infile   = os.path.join(temp, filename)
            outfile  = os.path.join(self._destination, filename)
        
            fileutils.move_file(infile, outfile)
