#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import copy

import pysam

import paleomix.common.fileutils as fileutils
import paleomix.common.utilities as utilities

from paleomix.common.formats.fasta import FASTA
from paleomix.common.formats.msa import MSA
from paleomix.node import NodeError, Node


class CollectSequencesNode(Node):
    def __init__(self, fasta_files, sequences, destination, dependencies=()):
        """
        fasta_files -- { taxon_name_1 : filename_1, ... }
        sequences   -- { interval_name_1, ... }
        """

        self._infiles = copy.deepcopy(fasta_files)
        self._sequences = utilities.safe_coerce_to_frozenset(sequences)
        self._destination = copy.copy(destination)
        self._outfiles = [
            os.path.join(destination, name + ".fasta") for name in self._sequences
        ]

        input_files = list(self._infiles.values())
        for filename in self._infiles.values():
            input_files.append(filename + ".fai")

        desc = "<CollectSequences: %i sequences from %i files -> '%s'>" % (
            len(self._sequences),
            len(self._infiles),
            self._destination,
        )
        Node.__init__(
            self,
            description=desc,
            input_files=input_files,
            output_files=self._outfiles,
            dependencies=dependencies,
        )

    def _setup(self, _config, _temp):
        for filename in self._infiles.values():
            with open(filename + ".fai") as handle:
                sequences = set()
                for line in handle:
                    sequences.add(line.split("\t", 1)[0])

                missing_sequences = list(self._sequences - sequences)
                if missing_sequences:
                    if len(missing_sequences) >= 4:
                        missing_sequences = missing_sequences[:3]
                        missing_sequences.append("...")

                    message = (
                        "FASTA file does not contain expected "
                        "sequences:\n  File =  %r\n  "
                        "Sequences = %s\n"
                    ) % (filename, ", ".join(missing_sequences))
                    raise NodeError(message)

    def _run(self, _config, temp):
        fasta_files = []
        for (name, filename) in sorted(self._infiles.items()):
            fasta_files.append((name, pysam.FastaFile(filename)))

        for sequence_name in sorted(self._sequences):
            filename = os.path.join(temp, sequence_name + ".fasta")
            with open(filename, "w") as out_handle:
                for (sample, fasta_file) in fasta_files:
                    sequence = fasta_file.fetch(sequence_name)
                    fasta = FASTA(sample, sequence_name, sequence)
                    fasta.write(out_handle)

    def _teardown(self, _config, temp):
        for destination in sorted(self._outfiles):
            source = fileutils.reroot_path(temp, destination)
            fileutils.move_file(source, destination)


class FilterSingletonsNode(Node):
    def __init__(self, input_file, output_file, filter_by, dependencies):
        self._input_file = input_file
        self._output_file = output_file
        self._filter_by = dict(filter_by)
        for (to_filter, groups) in self._filter_by.items():
            # The taxa to be filtered is implied to be part of the group,
            # but is not needed when actually carrying out the filtering
            groups = utilities.safe_coerce_to_frozenset(
                groups
            ) - utilities.safe_coerce_to_frozenset(to_filter)

            if not groups:
                raise RuntimeError(
                    "Singleton filtering must involve at least one other taxa"
                )
            self._filter_by[to_filter] = groups

        Node.__init__(
            self,
            description="filtering singleton nucleotides in %s" % (input_file,),
            input_files=[input_file],
            output_files=[output_file],
            dependencies=dependencies,
        )

    def _run(self, _config, temp):
        alignment = MSA.from_file(self._input_file)
        for (to_filter, groups) in self._filter_by.items():
            alignment = alignment.filter_singletons(to_filter, groups)

        temp_filename = fileutils.reroot_path(temp, self._output_file)
        with open(temp_filename, "w") as handle:
            alignment.to_file(handle)
        fileutils.move_file(temp_filename, self._output_file)
