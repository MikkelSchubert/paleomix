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
import copy
import itertools
import collections

import pypeline.common.fileutils as fileutils
import pypeline.common.utilities as utilities

from pypeline.common.formats.fasta import \
     FASTA
from pypeline.common.formats.msa import \
     MSA
from pypeline.node import \
     NodeError, \
     Node, \
     MetaNode



class CollectSequencesNode(Node):
    def __init__(self, fasta_files, sequences, destination, dependencies = ()):
        """
        fasta_files -- { taxon_name_1 : filename_1, ... }
        sequences   -- { interval_name_1, ... }
        """

        self._infiles     = copy.deepcopy(fasta_files)
        self._sequences   = utilities.safe_coerce_to_frozenset(sequences)
        self._destination = copy.copy(destination)
        self._outfiles    = [os.path.join(destination, name + ".fasta") for name in self._sequences]

        Node.__init__(self,
                      description  = "<CollectSequences: %i sequences from %i files -> '%s'>" \
                            % (len(self._sequences), len(self._infiles), self._destination),
                      input_files  = self._infiles.values(),
                      output_files = self._outfiles,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        fasta_names = []
        fasta_files = []
        for (name, filename) in self._infiles.iteritems():
            fasta_names.append(name)
            fasta_files.append(FASTA.from_file(filename))

        # Names of sequences that have been processed / are left to processes
        seq_names_done  = set()
        seq_names_left  = set(self._sequences)
        expected_num_records = len(self._infiles)
        def _process_records(sequence_name, records):
            if len(records) != expected_num_records:
                raise NodeError("Did not find expected number of records for %r; Expected %i, found %i" \
                                % (sequence_name, expected_num_records, len(records)))
            elif sequence_name in seq_names_done:
                raise NodeError("Multiple sequences with the same name: %r"\
                                % (sequence_name,))
            elif sequence_name in seq_names_left:
                seq_names_done.add(sequence_name)
                seq_names_left.remove(sequence_name)
                self._write_fasta(temp, fasta_names, sequence_name, records)

        fasta_partial = collections.defaultdict(dict)
        for records in itertools.izip(*fasta_files):
            sequence_name = records[0].name
            # Handle unsorted sets of files
            if not all((record.name == sequence_name) for record in records):
                for (name, record) in zip(fasta_names, records):
                    if (record.name in seq_names_done) or (name in fasta_partial.get(record.name, ())):
                        raise NodeError("Multiple sequences with the same name: %r"\
                                        % (record.name,))
                    elif record.name in seq_names_left:
                        fasta_partial[record.name][name] = record
                continue

            _process_records(sequence_name, records)

        for (sequence_name, records) in fasta_partial.iteritems():
            records = [records[name] for name in fasta_names]
            _process_records(sequence_name, records)


    def _teardown(self, _config, temp):
        for sequence in self._sequences:
            filename = sequence + ".fasta"
            infile   = os.path.join(temp, filename)
            outfile  = os.path.join(self._destination, filename)

            fileutils.move_file(infile, outfile)


    @classmethod
    def _write_fasta(cls, temp, fasta_names, sequence_name, records):
        filename = os.path.join(temp, sequence_name + ".fasta")
        with open(filename, "w") as out_handle:
            for (name, record) in zip(fasta_names, records):
                fasta = FASTA(name, record.name, record.sequence)
                out_handle.write(str(fasta))




class FilterSingletonsNode(Node):
    def __init__(self, input_file, output_file, filter_by, dependencies):
        self._input_file      = input_file
        self._output_file     = output_file
        self._filter_by       = dict(filter_by)
        for (to_filter, groups) in self._filter_by.items():
            if not groups:
                raise RuntimeError("Singleton filtering must involve at least one taxa")
            self._filter_by[to_filter] = groups

        Node.__init__(self,
                      description  = "<FilterSingleton: '%s' -> '%s'>" \
                            % (input_file, output_file),
                      input_files  = [input_file],
                      output_files = [output_file],
                      dependencies = dependencies)

    def _run(self, _config, temp):
        alignment = MSA.from_file(self._input_file)
        for (to_filter, groups) in self._filter_by.iteritems():
            alignment = alignment.filter_singletons(to_filter, groups)

        temp_filename = fileutils.reroot_path(temp, self._output_file)
        with open(temp_filename, "w") as handle:
            alignment.to_file(handle)
        fileutils.move_file(temp_filename, self._output_file)




class FilterSingletonsMetaNode(MetaNode):
    def __init__(self, input_files, destination, filter_by, dependencies = ()):
        subnodes = []
        filter_by = dict(filter_by)
        for (filename, node) in input_files.iteritems():
            output_filename = fileutils.reroot_path(destination, filename)
            subnodes.append(FilterSingletonsNode(input_file   = filename,
                                                 output_file  = output_filename,
                                                 filter_by    = filter_by,
                                                 dependencies = node))

        MetaNode.__init__(self,
                          description  = "<FilterSingleton: %i files -> '%s'>" \
                            % (len(subnodes), destination),
                          subnodes     = subnodes,
                          dependencies = dependencies)
