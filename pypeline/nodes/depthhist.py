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
import sys
import datetime
import subprocess

import pysam

from pypeline.node import Node
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.common.text import padded_table
from pypeline.common.fileutils import reroot_path, move_file, swap_ext, describe_files
from pypeline.common.utilities import get_in, set_in, safe_coerce_to_tuple
from pypeline.nodes.picard import concatenate_input_bams
from pypeline.common.timer import BAMTimer


_MAX_DEPTH = 200
_MAX_CONTIGS = 100


class DepthHistogramNode(Node):
    def __init__(self, config, target_name, input_files, output_file, intervals_file = None, print_stats = False, max_contigs = _MAX_CONTIGS, dependencies = ()):
        self._target_name = target_name
        self._input_files = safe_coerce_to_tuple(input_files)
        self._output_file = output_file
        self._intervals   = intervals_file
        self._print_stats = print_stats
        self._max_contigs = max_contigs
        self._max_contigs_reached = False

        input_files = []
        input_files.extend(self._input_files)
        input_files.extend(swap_ext(input_file, ".bai") for input_file in self._input_files)
        if intervals_file:
            input_files.append(intervals_file)

        executables = ["coverageBed"] if intervals_file else ["genomeCoverageBed"]
        auxiliary_files = []
        for cmd in concatenate_input_bams(config, self._input_files)[0]:
            executables.extend(cmd.executables)
            auxiliary_files.extend(cmd.auxiliary_files)

        Node.__init__(self,
                      description  = "<DepthHistogram: %s -> '%s'>" \
                        % (describe_files(self._input_files),
                           self._output_file),
                      input_files  = input_files,
                      output_files = self._output_file,
                      dependencies = dependencies,
                      executables  = executables,
                      auxiliary_files = auxiliary_files)


    def _setup(self, config, temp):
        pipe = os.path.join(temp, "input_file")
        # Dictionary of named pipes created by the node
        self._pipes = {"input_file" : pipe}
        # Dictionary of lists of processes created by the node
        self._procs = {}
        # Dictionary of output tables created by the node
        self._tables = {}
        # Dictionary of (filename, file-handle) for output tables
        self._handle = {}

        if len(self._input_files) == 1:
            filename = iter(self._input_files).next()
            os.symlink(os.path.abspath(filename), pipe)
        else:
            os.mkfifo(pipe)
            procs, _ = concatenate_input_bams(config, self._input_files, out = pipe)
            for proc in procs:
                proc.run(temp)
            self._procs["cat"] = procs


    def _run(self, config, temp):
        region_names = self._create_tables(config, temp)

        table = {}
        for (key, (filename, handle)) in self._tables.iteritems():
            handle.close()
            self._read_table(key, table, filename)

        temp_filename = reroot_path(temp, self._output_file)
        self._write_table(table, temp_filename, region_names)


    def _teardown(self, _config, temp):
        temp_filename = reroot_path(temp, self._output_file)
        move_file(temp_filename, self._output_file)

        for filename in self._pipes.itervalues():
            os.remove(filename)

        for (filename, _) in self._tables.itervalues():
            os.remove(filename)

        intervals = os.path.join(temp, "intervals.bed")
        if os.path.exists(intervals):
            os.remove(intervals)

        for proc in self._procs.get("cat", ()):
            proc.commit(temp)

        if not self._print_stats:
            os.remove(os.path.join(temp, "pipe_coverage_%i.stdout" % id(self)))


    def _create_tables(self, _config, temp):
        # Opening pipe/symlink created in _setup()
        out = sys.stderr
        if not self._print_stats:
            out = open(os.path.join(temp, "pipe_coverage_%i.stdout" % id(self)), "w")

        with pysam.Samfile(self._pipes["input_file"]) as samfile:
            timer = BAMTimer(samfile, out = out)
            intervals, region_names = self._get_intervals(temp, samfile)
            mapping   = self._open_handles(temp, samfile, intervals)
            for (counter, read) in enumerate(samfile):
                if counter % 100000 == 0:
                    if self._processes_terminated():
                        break

                # 0xf04 = 0x800 | 0x400 | 0x200 | 0x100 | 0x004
                #         0x800 = Secondary alignment
                #         0x400 = PCR Duplicate
                #         0x200 = Failed QC
                #         0x100 = Alternative alignment
                #         0x004 = Unmapped read
                if read.flag & 0xf04:
                    continue

                try:
                    readgroup = read.opt("RG")
                except KeyError:
                    readgroup = None

                for handle in mapping[readgroup]:
                    handle.write(read)
                timer.increment(read = read)
            timer.finalize()

        if not self._print_stats:
            out.close()

        for handle in self._handle.itervalues():
            handle.close()

        for proclst in self._procs.itervalues():
            for proc in proclst:
                if proc.wait() != 0:
                    raise RuntimeError("Error while running process: Return-code = %i" % proc.wait())

        return region_names


    def _calc_cumfrac(self, counts):
        cumfrac = [0] * len(counts)
        running_total = sum(counts)
        total = float(running_total or 1)
        for (index, count) in enumerate(counts):
            cumfrac[index] = "%.4f" % (running_total / total,)
            running_total -= count

        return cumfrac[1:]


    def _calc_max_depth(self, counts):
        counts = counts[1:]
        running_total = sum(counts)
        if not running_total:
            return "NA"

        total = float(running_total)
        for (index, count) in enumerate(counts):
            if running_total / total < 0.005:
                # Index corresponds to the previous depth
                return index
            running_total -= count

        return "NA"


    def _write_table(self, table, filename, region_names):
        if not any(table["<NA>"][None][None][1:]):
            table.pop("<NA>")

        zeros = [0] * (_MAX_DEPTH + 1)
        first = lambda pair: (pair[0] or "")
        rows = [["Name", "Sample", "Library", "Contig", "Size", "MaxDepth"] + ["MD_%03i" % i for i in xrange(1, _MAX_DEPTH + 1)]]
        for sample, libraries in sorted(table.iteritems(), key = first):
            for library, regions in sorted(libraries.iteritems(), key = first):
                for (region, size) in sorted(region_names.items(), key = first):
                    key = (self._target_name, sample, library, region)
                    row = [("*" if value is None else value) for value in key]

                    counts = regions.get(region, zeros)
                    row.append(size)
                    row.append(self._calc_max_depth(counts))
                    row.extend(self._calc_cumfrac(counts))
                    rows.append(row)

                    if (key[-1] is None) and self._max_contigs_reached:
                        row = list(row)
                        row[3] = "<Genome>"
                        rows.append(row)

                rows.append("#")
            rows.extend("#")

        while rows and rows[-1] == "#":
            rows.pop()

        with open(filename, "w") as table_file:
            table_file.write(_HEADER % datetime.datetime.now().isoformat())
            for line in padded_table(rows):
                table_file.write(line)
                table_file.write("\n")


    def _read_table(self, key, table, filename):
        all_key = "all" if self._intervals else "genome"
        with open(filename) as handle:
            for line in handle:
                fields = line.split("\t")
                # 'all' is generated by coverageBed, as a catchall group
                if (fields[0] == all_key) and (len(fields) == 5):
                    name = None
                elif self._max_contigs_reached:
                    continue
                elif (len(fields) > (3 + 4)):
                    # Probably a BED6 file, 4 columns are from bedTools
                    name = fields[3]
                else:
                    name = fields[0]

                ckey = key + (name,)
                if not get_in(table, ckey, None):
                    set_in(table, ckey, [0] * (_MAX_DEPTH + 1))

                depth = min(_MAX_DEPTH, int(fields[-4]))
                get_in(table, ckey)[depth] += int(fields[-3])


    def _processes_terminated(self):
        for proclist in self._procs.itervalues():
            for proc in proclist:
                if isinstance(proc, AtomicCmd):
                    if proc.ready():
                        return True
                elif proc.poll():
                    return True
        return False


    def _open_handles(self, temp, samfile, intervals):
        mapping    = {}
        readgroups = self._get_readgroups(samfile)
        for (rg_id, rg) in readgroups.iteritems():
            mapping[rg_id] = set([self._create_handle(temp, samfile, intervals, None,     None),
                                  self._create_handle(temp, samfile, intervals, rg["SM"], None),
                                  self._create_handle(temp, samfile, intervals, rg["SM"], rg["LB"])])

        return mapping


    def _create_handle(self, temp, samfile, intervals, *key):
        key = tuple(key)
        if key not in self._handle:
            input_file  = os.path.join(temp, "input_pipe_%02i" % len(self._pipes))
            output_file = os.path.join(temp, "output_table_%02i.txt" % len(self._handle))
            os.mkfifo(input_file)

            self._pipes[key]  = input_file
            self._tables[key] = (output_file, open(output_file, "w"))
            if self._intervals:
                call = ("coverageBed", "-hist", "-abam", input_file, "-b", intervals)
            else:
                call = ("genomeCoverageBed", "-max", str(_MAX_DEPTH), "-ibam", input_file, "-g", intervals)

            self._procs[key]  = [subprocess.Popen(call,
                                                  stdout    = self._tables[key][-1],
                                                  preexec_fn = os.setsid,
                                                  close_fds = True)]

            self._handle[key] = pysam.Samfile(input_file, "wbu", template = samfile)
        return self._handle[key]


    def _get_intervals(self, temp, samfile):
        region_names = {}
        intervals = os.path.join(temp, "intervals.bed")
        with open(intervals, "w") as handle:
            if self._intervals:
                with open(self._intervals) as source:
                    for line in source:
                        fields = line.rstrip("\r\n").split("\t")[:6]
                        if len(fields) < 4:
                            assert len(fields) == 3
                            fields.extend((fields[0] + "*", "0", "+"))
                        handle.write("\t".join(fields) + "\n")

                        name   = fields[3]
                        length = int(fields[2]) - int(fields[1])
                        region_names[name] = region_names.get(name, 0) + length
            else:
                for seq in samfile.header["SQ"]:
                    handle.write("{SN}\t{LN}\n".format(**seq))
                    region_names[seq["SN"]] = seq["LN"]

                if len(region_names) > self._max_contigs:
                    self._max_contigs_reached = True
                    region_names = {}
            region_names[None] = sum(region_names.itervalues())

            return intervals, region_names


    @classmethod
    def _get_readgroups(cls, bamfile):
        readgroups = {None : {"ID" : None, "SM" : "<NA>", "LB" : "<NA>"}}
        for readgroup in bamfile.header.get('RG', []):
            readgroups[readgroup["ID"]] = readgroup

        return readgroups



_HEADER = \
"""# Timestamp: %s
#
# Columns:
#   Contig:   Contig, chromosome, or feature for which a depth histogram was
#             created. Unnamed features are named after the chromosome or
#             contig on which they are located, with a star appended. For
#             example "chr1*".
#   Size:     The total size of the region. Multiple features with the same
#             name are combined into one row, with the size representing to
#             total of these. Note that overlapping bases are counted 2 (or
#             more) times.
#   MaxDepth: Maximum depth to use when calling SNPs, in order to exclude
#             (at least) the 0.5%% most extreme sites based on read depth,
#             not including sites with depth 0.
#   MD_*:     Fraction of sites with a minimum depth of 1-200.
#
"""
