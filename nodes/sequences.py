from __future__ import with_statement

import textwrap

import pysam

from pypeline.node import Node


class CollectSequencesNode(Node):
    def __init__(self, infiles, outfile, sequence, dependencies):
        self._infiles = dict(infiles)
        self._outfile = outfile
        self._sequence = sequence

        Node.__init__(self,
                      description  = "<CollectSequences: '%s' from %i files -> '%s'>" \
                            % (sequence, len(infiles), outfile),
                      input_files  = self._infiles.values(),
                      output_files = [self._output_file],
                      dependencies = dependencies)


    def _run(self, config, temp):
        filename = os.path.join(temp, os.path.basename(self._outfile))

        with open(filename, "w") as fasta:
            for (name, source) in self._infiles.iteritems():
                handle = pysam.Fastafile(source)
                sequence = textwrap.fill(handle.fetch(self._sequence), 60)
                handle.close()
                
                fasta.write(">%s\n%s\n" % (name, sequence))

        return True


    def _teardown(self, config, temp):
        filename = os.path.join(temp, os.path.basename(self._outfile))
        os.rename(filename, self._outfile)

        return True
