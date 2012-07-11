#!/usr/bin/python

import os
import itertools
import collections

from pypeline.node import Node, NodeError
from pypeline.fileutils import move_file, reroot_path

from pypeline.pylib.sequences import read_fasta
from pypeline.pylib.utilities import grouper, safe_coerce_to_tuple



class FastaToPhyNode(Node):
    def __init__(self, infiles, out_phy, out_partitions = None, strict = False, add_option_flag = False, dependencies = ()):
        self._infiles   = safe_coerce_to_tuple(infiles)
        self._add_flag  = add_option_flag
        self._strict    = strict
        self._out_phy   = out_phy
        self._out_part  = out_partitions

        description  = "<FastaToPhy: %s -> '%s'%s%s%s>" % \
            (("%i files" % len(infiles)) if (len(infiles) > 1) else ("'%s'" % infiles[0]),
             out_phy,
             ", strict" if strict else "",
             ", partitions" if out_partitions else "",
             ", options flag" if add_option_flag else "")
            
        output_files = [out_phy]
        if out_partitions:
            output_files.append(out_partitions)

        Node.__init__(self, 
                      description  = description,
                      input_files  = infiles,
                      output_files = output_files,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        by_group, by_sequence = self._read_sequences(self._infiles)

        # Length that can be handled by PHYLIP, PAML, etc.
        name_length = 10
        if not self._strict:
            # Length that can be handled by PhyML, RAxML, etc.
            name_length = min(30, max(len(group) for group in by_group))
        
        _name_format = '{0:<%i}' % name_length
        def format_name(group):
            return _name_format.format(group)[:name_length]

        self._write_sequences_interleaved(by_group, format_name, temp)
        if self._out_part:
            self._write_partitions(by_sequence, temp)


    def  _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._out_phy), self._out_phy)
        if self._out_part:
            move_file(reroot_path(temp, self._out_part), self._out_part)


    def _write_sequences_interleaved(self, by_group, format_name, temp):
        streams = []
        for (group, sequences) in sorted(by_group.iteritems()):
            sequence   = "".join((sequence for (_, sequence) in sorted(sequences.iteritems())))
            sequence_length = len(sequence)
            
            name       = format_name(group)
            row        = name
            row       +=  ("  " if ((len(name) - 9) % 12) else " ")
            
            frag_size  = min(10, max(0, 10 - (len(name) + 2) % 12))
            row       += sequence[:frag_size]
            sequence   = sequence[frag_size:]
            
            while (len(row) < 70) and sequence:
                row += "  " + sequence[:10]
                sequence = sequence[10:]

            rows = [row]
            for fragment in grouper(60, sequence, fillvalue = ""):
                rows.append(("  ".join("".join(column) for column in grouper(10, fragment, fillvalue = ""))))

            streams.append(rows)
        
        with open(reroot_path(temp, self._out_phy), "w") as output:
            output.write("   %i %i%s\n" % (len(by_group), sequence_length, " I" if self._add_flag else ""))
            for lst in itertools.izip(*streams):
                for row in lst:
                    output.write(row + "\n")
                output.write("\n")


    def _write_partitions(self, by_sequence, temp):
        end_pos = 0
        with open(reroot_path(temp, self._out_part), "w") as partitions:
            for (sequence_name, groups) in sorted(by_sequence.iteritems()):
                length = len(groups.itervalues().next())

                start_pos, end_pos = end_pos + 1, end_pos + length
        
                partitions.write("DNA, %s_12 = %i-%i\\3, %i-%i\\3\n" \
                                     % (sequence_name, start_pos, end_pos, start_pos + 1, end_pos))
                partitions.write("DNA, %s_3 = %i-%i\\3\n" \
                                     % (sequence_name, start_pos + 2, end_pos))                


    @classmethod
    def _read_sequences(cls, filenames):
        by_group    = collections.defaultdict(dict)
        by_sequence = collections.defaultdict(dict)

        for filename in sorted(filenames):
            sequence_name = os.path.splitext(os.path.basename(filename))[0]

            for (group, sequence) in read_fasta(filename):
                if (group in by_sequence[sequence_name]) or (sequence_name in by_group[group]):
                    raise NodeError("Multiple sequences with the same name (%s) for group '%s'." \
                                        % (sequence_name, group))

                by_group[group][sequence_name]    = sequence
                by_sequence[sequence_name][group] = sequence
                
        return cls._validate_sequences(by_group, by_sequence)


    @classmethod
    def _validate_sequences(cls, by_group, by_sequence):
        expected_groups = set(by_group)
        for (sequence_name, groups) in by_sequence.iteritems():
            missing_or_extra_groups = expected_groups.symmetric_difference(groups)
            if missing_or_extra_groups:
                raise NodeError("Unexpected/missing groups for sequence (%s): %s" \
                                    % (sequence_name, ", ".join(missing_or_extra_groups)))
        
            sequence_lengths = set(len(sequences) for sequences in groups.itervalues())
            if len(sequence_lengths) > 1:
                sequence_lengths = [str(length) for length in sequence_lengths]
                raise NodeError("Different lengths observed for the same sequence (%s): %s" \
                                    % (sequence_name, ", ".join(sequence_lengths)))

        return by_group, by_sequence
