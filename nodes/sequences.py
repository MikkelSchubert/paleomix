from __future__ import with_statement

import os
import textwrap

import pysam

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
        if not os.path.exists(self._destination):
            os.makedirs(self._destination)

        for sequence in self._sequences:
            filename = sequence + ".fasta"
            infile   = os.path.join(temp, filename)
            outfile  = os.path.join(self._destination, filename)
        
            os.rename(infile, outfile)
