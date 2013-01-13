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
from pypeline.common.text import padded_table, parse_padded_table
from pypeline.common.samwrap import SamfileReader
from pypeline.common.fileutils import reroot_path, move_file, swap_ext
from pypeline.common.utilities import get_in, set_in, safe_coerce_to_tuple


class CoverageNode(Node):
    def __init__(self, input_file, name, collapsed_prefix = "M_", dependencies = ()):
        self._name = name
        self._output_file = swap_ext(input_file, ".coverage")
        self._collapsed_prefix = collapsed_prefix

        Node.__init__(self,
                      description  = "<Coverage: '%s' -> '%s'>" \
                          % (input_file, self._output_file),
                      input_files  = input_file,
                      output_files = self._output_file,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        with SamfileReader(self.input_files[0]) as bamfile:
            raw_tables = self.initialize_raw_tables(bamfile)
            self.read_records(bamfile, raw_tables, self._collapsed_prefix)
            table = self.collapse_raw_tables(bamfile, raw_tables, self._name)

        _write_table(table, reroot_path(temp, self._output_file))
        move_file(reroot_path(temp, self._output_file), self._output_file)


    @classmethod
    def initialize_raw_tables(cls, bamfile):
        tables = {}
        for rg in bamfile.header['RG']:
            subtables = []
            for _ in bamfile.lengths:
                subtables.append({"SE" : 0, "PE_1" : 0, "PE_2" : 0, "Collapsed" : 0, "Hits" : 0, "M" : 0, "I" : 0, "D" : 0})
            tables[rg["ID"]] = subtables

        return tables


    @classmethod
    def read_records(cls, bamfile, tables, collapsed_prefix):
        for record in bamfile:
            if record.is_unmapped:
                continue

            readgroup = dict(record.tags)["RG"]
            subtable  = tables[readgroup][record.tid]

            if collapsed_prefix and record.qname.startswith(collapsed_prefix):
                subtable["Collapsed"] += 1
            else:
                flag = record.flag
                if flag & 0x40: # first of pair
                    subtable["PE_1"] += 1
                elif flag & 0x80: # second of pair
                    subtable["PE_2"] += 1
                else: # Singleton
                    subtable["SE"] += 1

            for (op, num) in record.cigar:
                if op < 3:
                    subtable["MID"[op]] += num

    @classmethod
    def collapse_raw_tables(cls, bamfile, tables, name):
        rg_to_library = {}
        rg_to_sample = {}
        for rg in bamfile.header['RG']:
            rg_to_library[rg["ID"]] = rg["LB"]
            rg_to_sample[rg["ID"]] = rg["SM"]

        collapsed = {}
        for rg in bamfile.header['RG']:
            subtables = []
            for size in bamfile.lengths:
                subtables.append({"SE" : 0, "PE_1" : 0, "PE_2" : 0, "Collapsed" : 0, "Size" : size, "Hits" : 0, "M" : 0, "I" : 0, "D" : 0})
            set_in(collapsed, (name, rg["SM"], rg["LB"]), subtables)

        for (readgroup, contigs) in tables.iteritems():
            library      = rg_to_library[readgroup]
            sample       = rg_to_sample[readgroup]
            subtable     = get_in(collapsed, (name, sample, library))
            for (contig_id, subtable_src) in enumerate(contigs):
                subtable_dst = subtable[contig_id]
                for (key, value) in subtable_src.iteritems():
                    subtable_dst[key] += value

        for (name, samples) in collapsed.items():
            for (sample, libraries) in samples.items():
                for (library, contigs) in libraries.items():
                    collapsed[name][sample][library] = dict(zip(bamfile.references, contigs))

        return collapsed



class MergeCoverageNode(Node):
    def __init__(self, input_files, output_file, dependencies = ()):
        self._input_files = safe_coerce_to_tuple(input_files)
        self._output_file = output_file

        Node.__init__(self,
                      description  = "<MergeCoverage: '%s' -> '%s'>" \
                          % (self._desc_files(self._input_files), self._output_file),
                      input_files  = self._input_files,
                      output_files = self._output_file,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        table = {}
        for filename in self._input_files:
            read_table(table, filename)

        _write_table(table, reroot_path(temp, self._output_file))
        move_file(reroot_path(temp, self._output_file), self._output_file)



def _calculate_totals(table):
    for (name, samples) in sorted(table.items()):
        set_in(table, (name, "*", "*", "*"), _calculate_totals_in(table))
        for (sample, libraries) in sorted(samples.items()):
            set_in(table, (name, sample, "*", "*"), _calculate_totals_in(libraries))
            for (library, contigs) in sorted(libraries.items()):
                set_in(table, (name, sample, library, "*"), _calculate_totals_in(contigs))
            set_in(table, (name, sample, "*", "*", "Size"), get_in(table, (name, sample, library, "*", "Size")))
        set_in(table, (name, "*", "*", "*", "Size"), get_in(table, (name, sample, "*", "*", "Size")))



def _calculate_totals_in(tables):
    totals = {"SE" : 0, "PE_1" : 0, "PE_2" : 0, "Collapsed" : 0, "Size" : 0, "Hits" : 0, "M" : 0, "I" : 0, "D" : 0}
    subtables = tables.values()
    while subtables:
        for subtable in list(subtables):
            if "SE" in subtable:
                for key in subtable:
                    totals[key] += subtable[key]
            else:
                subtables.extend(subtable.values())
            subtables.remove(subtable)
    return totals


def _build_rows(table):
    rows = [("Name", "Sample", "Library", "Contig", "Size", "Hits", "SE", "PE_1", "PE_2", "Collapsed", "M", "I", "D", "Coverage")]
    for (name, samples) in sorted(table.items()):
        for (sample, libraries) in sorted(samples.items()):
            for (library, contigs) in sorted(libraries.items()):
                for (contig, subtable) in sorted(contigs.items()):
                    row = [name,
                           sample,
                           library,
                           contig,
                           subtable["Size"],
                           subtable["SE"] + subtable["PE_1"] + subtable["PE_2"] + subtable["Collapsed"],
                           subtable["SE"],
                           subtable["PE_1"],
                           subtable["PE_2"],
                           subtable["Collapsed"],
                           subtable["M"],
                           subtable["I"],
                           subtable["D"],
                           float(subtable["M"]) / subtable["Size"]]
                    rows.append(row)
                rows.append("#")
            rows.append("#")

    while rows[-1] == "#":
        rows.pop()
    return rows


def read_table(table, filename):
    with open(filename) as table_file:
        for record in parse_padded_table(table_file):
            key = (record["Name"], record["Sample"], record["Library"], record["Contig"])
            if "*" in key:
                continue

            subtable = get_in(table, key, {"Size" : int(record["Size"])})
            assert int(subtable["Size"]) == int(record["Size"])

            for field in ("Hits", "SE", "PE_1", "PE_2", "Collapsed", "M", "I", "D"):
                subtable[field] = subtable.get(field, 0) + int(record.get(field, 0))
            set_in(table, key, subtable)


def _write_table(table, filename):
    _calculate_totals(table)
    rows = _build_rows(table)
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
#  Collapsed:      Number of hits for PE pair collapsed into a single read
#  M, I, D:        Number of aligned (M), inserted (I) and deleted (D) bases relative to references
#  Coverage:       Estimated coverage based on aligned (M) bases
"""
