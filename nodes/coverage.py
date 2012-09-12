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
import datetime
import collections

from pypeline.node import Node
from pypeline.common.text import padded_table
from pypeline.common.samwrap import SamfileReader
from pypeline.common.fileutils import reroot_path, move_file, swap_ext
from pypeline.common.utilities import get_in, set_in


class CoverageNode(Node):
    def __init__(self, input_file, name, dependencies = ()):
        self._name = name
        self._input_file  = input_file
        self._output_file = swap_ext(input_file, ".coverage")

        Node.__init__(self, 
                      description  = "<Coverage: '%s' -> '%s'>" \
                          % (input_file, self._output_file),
                      input_files  = self._input_file,
                      output_files = self._output_file,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        table = dict()
        with SamfileReader(self._input_file) as bamfile:
            rg_to_library = {}
            rg_to_sample = {}
            for rg in bamfile.header['RG']:
                rg_to_library[rg["ID"]] = rg["LB"]
                rg_to_sample[rg["ID"]] = rg["SM"]
            
            for record in bamfile:
                if record.is_unmapped:
                    continue

                readgroup = dict(record.tags)["RG"]
                library   = rg_to_library[readgroup]
                sample    = rg_to_sample[readgroup]
                contig    = bamfile.references[record.tid]
                key       = (sample, library, contig)
                
                subtable  = get_in(table, key)
                if not subtable:
                    subtable = {"SE" : 0, "PE_1" : 0, "PE_2" : 0, "CIGAR" : [0] * 9}
                    set_in(table, key, subtable)

                flag = record.flag
                if flag & 0x40: # first of pair
                    subtable["PE_1"] += 1
                elif flag & 0x80: # second of pair
                    subtable["PE_2"] += 1
                else: # Singleton
                    subtable["SE"] += 1

                cigar = subtable["CIGAR"]
                for (op, num) in record.cigar:
                    cigar[op] += num

            lengths = dict(zip(bamfile.references, bamfile.lengths))
            lengths["*"] = sum(lengths.values())

        set_in(table, "***", self._calculate_totals(table))
        for (sample, libraries) in sorted(table.items()):
            set_in(table, (sample, "*", "*"), self._calculate_totals(libraries))
            for (library, contigs) in sorted(libraries.items()):
                set_in(table, (sample, library, "*"), self._calculate_totals(contigs))
            
             
        rows = self._build_rows(table, lengths, self._name)
        self._write_rows(rows, reroot_path(temp, self._output_file))
        move_file(reroot_path(temp, self._output_file), self._output_file)

    
    @classmethod
    def _calculate_totals(cls, tables):
        totals = {"SE" : 0, "PE_1" : 0, "PE_2" : 0, "CIGAR" : [0] * 9}
        subtables = tables.values()
        while subtables:
            for subtable in list(subtables):
                if "SE" in subtable:
                    totals["SE"] += subtable["SE"]
                    totals["PE_1"] += subtable["PE_1"]
                    totals["PE_2"] += subtable["PE_2"]
                
                    for (index, count) in enumerate(subtable["CIGAR"]):
                        totals["CIGAR"][index] += count
                else:
                    subtables.extend(subtable.values())
                subtables.remove(subtable)
        
        return totals


    def _build_rows(cls, table, lengths, name):
        rows = [("Name", "Sample", "Library", "Contig", "Size", "Hits", "SE", "PE_1", "PE_2", "M", "I", "D", "Coverage")]
        for (sample, libraries) in sorted(table.items()):
            for (library, contigs) in sorted(libraries.items()):
                for (contig, subtable) in sorted(contigs.items()):
                    # Sum aligned bases ('M'), matches ('='), and mismatches ('X')
                    aligned_nts = subtable["CIGAR"][0] + subtable["CIGAR"][7] + subtable["CIGAR"][8]
                    row = [name,
                           sample,
                           library, 
                           contig, 
                           lengths[contig],
                           subtable["SE"] + subtable["PE_1"] + subtable["PE_2"],
                           subtable["SE"],
                           subtable["PE_1"],
                           subtable["PE_2"],
                           aligned_nts,
                           subtable["CIGAR"][1],
                           subtable["CIGAR"][2],
                           float(aligned_nts) / lengths[contig]]
                    rows.append(row)
                rows.append("#")
            rows.append("#")

        while rows[-1] == "#":
            rows.pop()
        
        return rows


    def _write_rows(cls, rows, filename):
        with open(filename, "w") as table_file:
            table_file.write(_HEADER % datetime.datetime.now().isoformat())
            for line in padded_table(rows):
                table_file.write(line)
                table_file.write("\n")



_HEADER = \
"""# Timestamp: %s
#
# Columns:
#  Hits:           Sum of SE, PE_1, and PE_2 hits
#  SE, PE_1, PE_2: Number of Single Ended, and Pair Ended (mate 1 and 2) hits
#  M, I, D:        Number of aligned (M), inserted (I) and deleted (D) bases relative to references
#  Coverage:       Estimated coverage based on aligned (M) bases
"""
