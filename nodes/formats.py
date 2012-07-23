#!/usr/bin/python

import os
import itertools
import collections

from pypeline.node import Node, NodeError
from pypeline.fileutils import move_file, reroot_path

from pypeline.common.formats.msa import read_msa, join_msa
from pypeline.common.formats.phylip import interleaved_phy
from pypeline.common.utilities import grouper, safe_coerce_to_tuple








class FastaToPartitionsNode(Node):
    def __init__(self, infiles, out_partitions, partition_by = "123", dependencies = ()):
        assert len(partition_by) == 3
        self._out_part  = out_partitions
        self._part_by   = partition_by

        description  = "<FastaToPartitions: %i file(s) -> '%s'>" % \
            (len(infiles), out_partitions)
            
        Node.__init__(self, 
                      description  = description,
                      input_files  = infiles,
                      output_files = out_partitions,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        partition_by = {}
        for (group, offset) in zip(self._part_by, range(3)):
            partition_by.setdefault(group, []).append(offset)
        partition_by = list(sorted(partition_by.items()))

        end = 0
        partitions = []
        for (name, msa) in _read_sequences(self.input_files):
            length = len(msa.itervalues().next())
            start, end = end + 1, end + length

            for (group, offsets) in partition_by:
                if len(offsets) != 3:
                    parts = [("%i-%i\\3" % (start + offset, end)) for offset in offsets]
                else:
                    parts = ["%i-%i" % (start, end)]

                partition = "DNA, %s_%s = %s\n" \
                    % (name, group, ", ".join(parts))
                partitions.append(partition)

        with open(reroot_path(temp, self._out_part), "w") as part_file:
            part_file.writelines(partitions)


    def  _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._out_part), self._out_part)




class FastaToInterleavedPhyNode(Node):
    def __init__(self, infiles, out_phy, add_flag = False, dependencies = ()):
        self._add_flag  = add_flag
        self._out_phy   = out_phy

        description  = "<FastaToInterleavedPhy: %i file(s) -> '%s'%s>" % \
            (len(infiles), out_phy, (" (w/ flag)" if add_flag else ""))
            
        Node.__init__(self, 
                      description  = description,
                      input_files  = infiles,
                      output_files = [out_phy],
                      dependencies = dependencies)


    def _run(self, _config, temp):
        msa = join_msa(*(read_msa(filename) for filename in self.input_files))

        with open(reroot_path(temp, self._out_phy), "w") as output:
            output.write(interleaved_phy(msa, add_flag = self._add_flag))


    def  _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._out_phy), self._out_phy)





def _read_sequences(filenames):
    expected_groups = None
    for filename in sorted(filenames):
        name = os.path.splitext(os.path.basename(filename))[0]
        msa  = read_msa(filename)

        if not expected_groups:
            expected_groups = set(msa)
        elif set(msa) != expected_groups:
            difference = expected_groups.symmetric_difference(msa)
            raise NodeError("Unexpected/missing groups for sequence (%s): %s" \
                                % (name, ", ".join(difference)))
        
        yield (name, msa)
