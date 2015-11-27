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
from __future__ import with_statement

import os
import copy
import itertools
import collections

import pysam

import paleomix.common.fileutils as fileutils
import paleomix.common.utilities as utilities
import paleomix.common.sequences as sequtils
import paleomix.common.text as text

from paleomix.common.formats.fasta import \
    FASTA
from paleomix.common.formats.msa import \
    MSA
from paleomix.node import \
    NodeError, \
    Node
from paleomix.common.bedtools import \
    BEDRecord


class CollectSequencesNode(Node):
    def __init__(self, fasta_files, sequences, destination, dependencies=()):
        """
        fasta_files -- { taxon_name_1 : filename_1, ... }
        sequences   -- { interval_name_1, ... }
        """

        self._infiles = copy.deepcopy(fasta_files)
        self._sequences = utilities.safe_coerce_to_frozenset(sequences)
        self._destination = copy.copy(destination)
        self._outfiles = [os.path.join(destination, name + ".fasta")
                          for name in self._sequences]

        input_files = list(self._infiles.itervalues())
        for filename in self._infiles.itervalues():
            input_files.append(filename + ".fai")

        desc = "<CollectSequences: %i sequences from %i files -> '%s'>" \
               % (len(self._sequences), len(self._infiles), self._destination)
        Node.__init__(self,
                      description=desc,
                      input_files=input_files,
                      output_files=self._outfiles,
                      dependencies=dependencies)

    def _setup(self, _config, _temp):
        for filename in self._infiles.itervalues():
            with open(filename + ".fai") as handle:
                sequences = set()
                for line in handle:
                    sequences.add(line.split("\t", 1)[0])

                missing_sequences = list(self._sequences - sequences)
                if missing_sequences:
                    if len(missing_sequences) >= 4:
                        missing_sequences = missing_sequences[:3]
                        missing_sequences.append("...")

                    message = ("FASTA file does not contain expected "
                               "sequences:\n  File =  %r\n  "
                               "Sequences = %s\n") \
                        % (filename, ", ".join(missing_sequences))
                    raise NodeError(message)

    def _run(self, _config, temp):
        fasta_files = []
        for (name, filename) in sorted(self._infiles.iteritems()):
            fasta_files.append((name, pysam.Fastafile(filename)))

        for sequence_name in sorted(self._sequences):
            filename = os.path.join(temp, sequence_name + ".fasta")
            with open(filename, "w") as out_handle:
                for (sample, fasta_file) in fasta_files:
                    sequence = fasta_file.fetch(sequence_name)
                    fasta = FASTA(sample, sequence_name, sequence)
                    out_handle.write(str(fasta))

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
            groups = utilities.safe_coerce_to_frozenset(groups) \
                - utilities.safe_coerce_to_frozenset(to_filter)

            if not groups:
                raise RuntimeError("Singleton filtering must involve at least "
                                   "one other taxa")
            self._filter_by[to_filter] = groups

        Node.__init__(self,
                      description="<FilterSingleton: '%s' -> '%s'>"
                      % (input_file, output_file),
                      input_files=[input_file],
                      output_files=[output_file],
                      dependencies=dependencies)

    def _run(self, _config, temp):
        alignment = MSA.from_file(self._input_file)
        for (to_filter, groups) in self._filter_by.iteritems():
            alignment = alignment.filter_singletons(to_filter, groups)

        temp_filename = fileutils.reroot_path(temp, self._output_file)
        with open(temp_filename, "w") as handle:
            alignment.to_file(handle)
        fileutils.move_file(temp_filename, self._output_file)


class ExtractReferenceNode(Node):
    def __init__(self, reference, bedfile, outfile, dependencies=()):
        self._reference = reference
        self._bedfile = bedfile
        self._outfile = outfile

        description = "<ExtractReference: '%s' -> '%s'>" \
            % (reference, outfile)
        Node.__init__(self,
                      description=description,
                      input_files=[reference, bedfile],
                      output_files=[outfile],
                      dependencies=dependencies)

    def _run(self, _config, temp):
        def _by_name(bed):
            return bed.name

        fastafile = pysam.Fastafile(self._reference)
        seqs = collections.defaultdict(list)
        with open(self._bedfile) as bedfile:
            bedrecords = text.parse_lines_by_contig(bedfile, BEDRecord)
            for (contig, beds) in sorted(bedrecords.iteritems()):
                beds.sort(key=lambda bed: (bed.contig, bed.name, bed.start))

                for (gene, gene_beds) in itertools.groupby(beds, _by_name):
                    gene_beds = tuple(gene_beds)
                    sequence = self._collect_sequence(fastafile, gene_beds)
                    seqs[(contig, gene)] = sequence

        temp_file = os.path.join(temp, "sequences.fasta")
        with open(temp_file, "w") as out_file:
            for ((_, gene), sequence) in sorted(seqs.items()):
                FASTA(gene, None, sequence).write(out_file)

        fileutils.move_file(temp_file, self._outfile)

    @classmethod
    def _collect_sequence(cls, fastafile, beds):
        sequence = []
        for bed in beds:
            fragment = fastafile.fetch(bed.contig, bed.start, bed.end)
            if len(fragment) != (bed.end - bed.start):
                cls._report_failure(bed, fragment)

            sequence.append(fragment)
        sequence = "".join(sequence)

        if any((bed.strand == "-") for bed in beds):
            assert all((bed.strand == "-") for bed in beds)
            sequence = sequtils.reverse_complement(sequence)

        return sequence

    @classmethod
    def _report_failure(cls, bed, fragment):
        message = "Failed to extract region from " \
                  "reference sequence at %s:%i-%i; got " \
                  "%i bp, but expected %i bp." \
                  % (bed.contig, bed.start, bed.end,
                     len(fragment), (bed.end - bed.start))
        raise NodeError(message)
