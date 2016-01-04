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
import re
import sys
import math
import numbers
import collections

from paleomix.node import Node, NodeError
from paleomix.common.utilities import set_in, get_in
from paleomix.common.fileutils import move_file, reroot_path
from paleomix.tools.bam_stats.coverage import \
    read_table as read_coverage_table
from paleomix.common.bedtools import BEDRecord

import paleomix.common.text as text


_PE_READS = frozenset(("Paired", "Singleton",
                       "Collapsed", "CollapsedTruncated"))
_SE_READS = frozenset(("Single",))
_BAMS     = frozenset(())


class SummaryTableNode(Node):
    def __init__(self, config, makefile, target, cov_for_lanes, cov_for_libs, dependencies = ()):
        self._target        = target.name
        self._output_file   = os.path.join(config.destination, self._target + ".summary")
        self._prefixes      = makefile["Prefixes"]
        self._makefile      = makefile["Statistics"]

        self._in_raw_bams = cov_for_lanes
        self._in_lib_bams = cov_for_libs
        input_files = set()
        input_files.update(sum(map(list, self._in_raw_bams.values()), []))
        input_files.update(sum(map(list, self._in_lib_bams.values()), []))

        self._in_raw_read = collections.defaultdict(list)
        for prefix in target.prefixes:
            for sample in prefix.samples:
                for library in sample.libraries:
                    for lane in library.lanes:
                        if lane.reads:
                            if lane.reads.stats:
                                value = lane.reads.stats
                                input_files.add(value)
                            elif set(lane.reads.files) & _PE_READS:
                                value = _PE_READS
                            elif set(lane.reads.files) & _SE_READS:
                                value = _SE_READS
                            else:
                                assert False
                        else:
                            value = _BAMS
                        self._in_raw_read[(sample.name, library.name, lane.name)] = value

        Node.__init__(self,
                      description  = "<Summary: %s>" % self._output_file,
                      input_files  = filter(None, input_files),
                      output_files = [self._output_file],
                      dependencies = dependencies)


    def _run(self, config, temp):
        rois    = self._stat_areas_of_interest(self._prefixes)
        genomes = self._stat_prefixes(self._prefixes)
        with open(reroot_path(temp, self._output_file), "w") as table:
            table.write("# Command:\n")
            table.write("#     %s\n" % (" ".join(sys.argv)),)
            table.write("#\n")
            table.write("# Directory:\n")
            table.write("#     %s\n" % (os.getcwd()),)
            table.write("#\n")
            table.write("# Makefile:\n")
            table.write("#     Filename: %s\n" % (self._makefile["Filename"],))
            table.write("#     SHA1Sum:  %s\n" % (self._makefile["Hash"],))
            table.write("#     MTime:    %s\n" % (self._makefile["MTime"],))
            table.write("#\n")
            self._write_genomes(table, genomes)
            table.write("#\n")
            self._write_areas_of_interest(table, rois)
            table.write("#\n#\n")

            for roi in rois.itervalues():
                genomes[roi["Label"]] = {"Size" : roi["Size"]}
            self._write_tables(table, genomes)


    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._output_file), self._output_file)


    def _write_genomes(self, table, genomes):
        table.write("# Genomes:\n")
        rows = [["Name", "Label", "Contigs", "Size", "Prefix"]]
        for (_, prefix) in sorted(self._prefixes.items()):
            stats = genomes[prefix["Name"]]
            rows.append((prefix["Name"], prefix.get("Label", "-"), stats["NContigs"], stats["Size"], prefix["Path"]))

        for line in text.padded_table(rows):
            table.write("#     %s\n" % (line,))


    def _write_areas_of_interest(self, table, rois):
        table.write("# Regions Of Interest:\n")
        rows = [["Genome", "ROI", "Size", "NFeatures", "NIntervals", "Path"]]
        for (_, roi) in sorted(rois.items()):
            rows.append([roi[key] for key in ("Genome", "Name", "Size", "NFeatures", "NIntervals", "Path")])

        for line in text.padded_table(rows):
            table.write("#     %s\n" % (line,))


    def _write_tables(self, out, genomes):
        rows = [["Target", "Sample", "Library", "Measure", "Value", "# Description"]]
        for (target, samples) in sorted(self._read_tables(self._prefixes, genomes).iteritems()):
            for (sample, libraries) in sorted(samples.iteritems()):
                for (library, prefixes) in sorted(libraries.iteritems()):
                    ordered = [("reads", prefixes.pop("reads"))] if "reads" in prefixes else []
                    ordered.extend(sorted(prefixes.items()))

                    for (prefix, table) in ordered:
                        table.pop("hits_unique_nts(%s)" % prefix, None)

                        for (key, (value, comment)) in sorted(table.iteritems(), key = _measure_ordering):
                            if isinstance(value, numbers.Number) and math.isnan(value):
                                value = "NA"
                            rows.append((target, sample, library, key, value, comment))
                        rows.append("")
                    rows.append("")

        for line in text.padded_table(rows):
            out.write("%s\n" % line)


    def _read_tables(self, prefixes, genomes):
        table = {}
        self._read_reads_settings(table)
        self._read_raw_bam_stats(table)
        self._read_lib_bam_stats(table)

        for (target, samples) in table.items():
            merged_samples = {}
            for (sample, libraries) in samples.items():
                merged_libraries = {}
                for (library, subtables) in libraries.items():
                    for (tblname, subtable) in subtables.items():
                        merged_libraries[tblname] = self._merge_tables((merged_libraries.get(tblname, {}), subtable))
                        merged_samples[tblname] = self._merge_tables((merged_samples.get(tblname, {}), subtable))
                    libraries[library] = self._annotate_subtables(subtables, genomes)
                set_in(table, (target, sample, "*"), self._annotate_subtables(merged_libraries, genomes))
            set_in(table, (target, "*", "*"), self._annotate_subtables(merged_samples, genomes))

        return table


    @classmethod
    def _annotate_subtables(cls, subtables, genomes):
        if "mitochondrial" in subtables and "nuclear" in subtables:
            subtables["endogenous"] = cls._create_endogenous_subtable(subtables, genomes)

        for (tblname, subtable) in subtables.iteritems():
            if tblname == "reads":
                fractions = [("seq_trash_se",      "seq_reads_se",       "seq_trash_se_frac",   "# Fraction of SE reads trashed"),
                             ("seq_trash_pe_1",    "seq_reads_pairs",    "seq_trash_pe_1_frac", "# Fraction of PE mate 1 reads trashed"),
                             ("seq_trash_pe_2",    "seq_reads_pairs",    "seq_trash_pe_2_frac", "# Fraction of PE mate 2 reads trashed"),
                             ("seq_collapsed",     "seq_reads_pairs",    "seq_collapsed_frac",  "# Fraction of PE pairs collapsed into one read"),
                             ("seq_retained_nts",  "seq_retained_reads", "seq_retained_length", "# Average number of NTs in retained reads")]

                for (numerator, denominator, measure, comment) in fractions:
                    if (numerator in subtable) and (denominator in subtable):
                        value = float(subtable[numerator][0]) / subtable[denominator][0]
                        subtable[measure] = (value, comment)
            else:
                total_hits  = subtable["hits_raw(%s)" % tblname][0]
                total_nts   = subtable["hits_unique_nts(%s)" % tblname][0]
                total_uniq  = subtable["hits_unique(%s)" % tblname][0]
                total_reads = subtables.get("reads",{}).get("seq_retained_reads", (float("NAN"),))[0]

                subtable["hits_raw_frac(%s)" % tblname] = (total_hits / float(total_reads), "# Total number of hits vs. total number of reads retained")
                subtable["hits_unique_frac(%s)" % tblname] = (total_uniq / float(total_reads), "# Total number of unique hits vs. total number of reads retained")
                subtable["hits_clonality(%s)" % tblname] = (1 - total_uniq / (float(total_hits) or float("NaN")), "# Fraction of hits that were PCR duplicates")
                subtable["hits_length(%s)" % tblname] = (total_nts / (float(total_uniq) or float("NaN")), "# Average number of aligned bases per unique hit")
                subtable["hits_coverage(%s)" % tblname] = (total_nts / float(genomes[tblname]["Size"]), "# Estimated coverage from unique hits")

        return subtables


    @classmethod
    def _create_endogenous_subtable(self, subtables, genomes):
        nucl = subtables["nuclear"]
        mito = subtables["mitochondrial"]

        total_hits            = mito["hits_raw(mitochondrial)"][0]         + nucl["hits_raw(nuclear)"][0]
        total_hits_unique     = mito["hits_unique(mitochondrial)"][0]      + nucl["hits_unique(nuclear)"][0]
        total_hits_unique_nts = mito["hits_unique_nts(mitochondrial)"][0]  + nucl["hits_unique_nts(nuclear)"][0]

        ratio_hits, ratio_genome, ratio_genome_inv = "NA", "NA", "NA"
        if mito["hits_unique(mitochondrial)"][0]:
            ratio_nts    = float(nucl["hits_unique_nts(nuclear)"][0]) / mito["hits_unique_nts(mitochondrial)"][0]
            ratio_hits   = float(nucl["hits_unique(nuclear)"][0]) / mito["hits_unique(mitochondrial)"][0]
            ratio_genome = ratio_nts / ((float(genomes["nuclear"]["Size"]) * 2) / float(genomes["mitochondrial"]["Size"]))
            ratio_genome_inv = ratio_genome ** -1

        return {
            "hits_raw(endogenous)"         : (total_hits,                       "# Total number of hits against the nuclear and mitochondrial genome"),
            "hits_unique(endogenous)"      : (total_hits_unique,                "# Total number of unique reads (PCR duplicates removed)"),
            "hits_unique_nts(endogenous)"  : (total_hits_unique_nts,            None),
            "ratio_reads(nuc,mito)"        : (ratio_hits,                       "# Ratio of unique hits: Hits(nuc) / H(mito)"),
            "ratio_genome(nuc,mito)"       : (ratio_genome,                     "# Ratio of NTs of unique hits corrected by genome sizes: (NTs(nuc) / NTs(mito)) / ((2 * Size(nuc)) / Size(mito))"),
            "ratio_genome(mito,nuc)"       : (ratio_genome_inv,                 "# Ratio of NTs of unique hits corrected by genome sizes: (NTs(mito) / NTs(nuc)) / (Size(mito) / (2 * Size(nuc)))")
            }


    def _read_reads_settings(self, table):
        for ((sample, library, barcode), filename) in self._in_raw_read.iteritems():
            key = (self._target, sample, library, "reads", barcode)
            set_in(table, key, self._stat_read_settings(filename))

        for (target, samples) in table.iteritems():
            for (sample, libraries) in samples.iteritems():
                for (library, prefixes) in libraries.iteritems():
                    prefixes["reads"] = self._merge_tables(prefixes["reads"].values())

        return table


    def _read_raw_bam_stats(self, table):
        for ((genome, target, sample, library), filenames) in self._in_raw_bams.iteritems():
            key = (target, sample, library)
            hits, _ = self._read_coverage_tables(key, filenames)

            value = (hits, "# Total number of hits (prior to PCR duplicate filtering)")
            set_in(table, (target, sample, library, genome, "hits_raw(%s)" % genome), value)


    def _read_lib_bam_stats(self, table):
        for ((genome, target, sample, library), filenames) in self._in_lib_bams.iteritems():
            key = (target, sample, library)
            hits, nts = self._read_coverage_tables(key, filenames)

            value = (hits, "# Total number of hits (excluding any PCR duplicates)")
            set_in(table, (target, sample, library, genome, "hits_unique(%s)" % genome), value)
            set_in(table, (target, sample, library, genome, "hits_unique_nts(%s)" % genome), (nts, None))


    @classmethod
    def _read_coverage_tables(cls, key, filenames):
        hits = nts = 0
        for filename in filenames:
            subtable = {}
            read_coverage_table(subtable, filename)
            contigtables = get_in(subtable, key)

            if contigtables is None:
                raise NodeError("Error reading table %r; row not found:"
                                "\n   %s   ...\n\nIf files have been renamed "
                                "during the run, then please remove this file "
                                "in that it may be re-generated.\nHowever, "
                                "note that read-group tags in the BAM files "
                                "may not be correct!"
                                % (filename, "   ".join(key)))

            for contigtable in contigtables.itervalues():
                hits += contigtable["Hits"]
                nts += contigtable["M"]
        return hits, nts


    @classmethod
    def _merge_tables(cls, tables):
        merged = {}
        for table in tables:
            for (measure, (value, comment)) in table.iteritems():
                if not isinstance(value, numbers.Number):
                    other, _ = merged.get(measure, (value, None))
                    merged[measure] = (value if (value == other) else "*", comment)
                else:
                    other, _ = merged.get(measure, (0, None))
                    merged[measure] = (value + other, comment)
        return merged


    @classmethod
    def _stat_read_settings(cls, filename):
        if isinstance(filename, frozenset):
            if (filename == _SE_READS):
                return {
                    "lib_type"           : ("SE", "# SE, PE, or * (for both)"),
                    "seq_reads_se"       : (float("nan"),  "# Total number of single-ended reads"),
                    "seq_trash_se"       : (float("nan"),  "# Total number of trashed reads"),
                    "seq_retained_nts"   : (float("nan"),  "# Total number of NTs in retained reads"),
                    "seq_retained_reads" : (float("nan"),  "# Total number of retained reads"),
                    }
            elif (filename == _PE_READS):
                return {
                    "lib_type"            : ("PE", "# SE, PE, or * (for both)"),
                    "seq_reads_pairs"     : (float("nan"), "# Total number of reads"),
                    "seq_trash_pe_1"      : (float("nan"), "# Total number of reads"),
                    "seq_trash_pe_2"      : (float("nan"), "# Total number of reads"),
                    "seq_retained_nts"    : (float("nan"), "# Total number of NTs in retained reads"),
                    "seq_retained_reads"  : (float("nan"), "# Total number of retained reads"),
                    "seq_collapsed"       : (float("nan"), "# Total number of pairs collapsed into one read"),
                    }
            else:
                return {
                    "lib_type"           : ("*", "# SE, PE, or * (for both)"),
                    "seq_reads_se"       : (float("nan"),  "# Total number of single-ended reads"),
                    "seq_trash_se"       : (float("nan"),  "# Total number of trashed reads"),
                    "seq_reads_pairs"    : (float("nan"), "# Total number of reads"),
                    "seq_trash_pe_1"     : (float("nan"), "# Total number of reads"),
                    "seq_trash_pe_2"     : (float("nan"), "# Total number of reads"),
                    "seq_retained_nts"   : (float("nan"),  "# Total number of NTs in retained reads"),
                    "seq_retained_reads" : (float("nan"),  "# Total number of retained reads"),
                    "seq_collapsed"      : (float("nan"), "# Total number of pairs collapsed into one read"),
                    }

        with open(filename) as settings_file:
            settings = settings_file.read()
            def _re_search(regexp, default = None):
                match = re.search(regexp, settings)
                if not match:
                    if default is not None:
                        return default
                    raise KeyError("Could not find match with RegExp %s in file '%s'" \
                                   % (repr(regexp), filename))

                return int(match.groups()[0])

            if "Paired end mode" in settings or "paired-end reads" in settings:
                return {
                    "lib_type"            : ("PE", "# SE, PE, or * (for both)"),
                    "seq_reads_pairs"     : (_re_search("number of read pairs: ([0-9]+)"),    "# Total number of pairs"),
                    "seq_trash_pe_1"      : (_re_search("discarded mate 1 reads: ([0-9]+)"),  "# Total number of reads"),
                    "seq_trash_pe_2"      : (_re_search("discarded mate 2 reads: ([0-9]+)"),  "# Total number of reads"),
                    "seq_retained_nts"    : (_re_search("retained nucleotides: ([0-9]+)"),    "# Total number of NTs in retained reads"),
                    "seq_retained_reads"  : (_re_search("retained reads: ([0-9]+)"),          "# Total number of retained reads"),
                    "seq_collapsed"       : (_re_search("of (?:full-length )?collapsed pairs: ([0-9]+)", 0) + \
                                             _re_search("of truncated collapsed pairs: ([0-9]+)", 0),
                                            "# Total number of pairs collapsed into one read"),
                    }
            elif "Single end mode" in settings or "single-end reads" in settings:
                return {
                    "lib_type"            : ("SE", "# SE, PE, or * (for both)"),
                    "seq_reads_se"        : (_re_search("number of (?:reads|read pairs): ([0-9]+)"),  "# Total number of single-ended reads"),
                    "seq_trash_se"        : (_re_search("discarded mate 1 reads: ([0-9]+)"),          "# Total number of trashed reads"),
                    "seq_retained_nts"    : (_re_search("retained nucleotides: ([0-9]+)"),            "# Total number of NTs in retained reads"),
                    "seq_retained_reads"  : (_re_search("retained reads: ([0-9]+)"),                  "# Total number of retained reads"),
                    }
            else:
                assert False, filename


    @classmethod
    def _stat_areas_of_interest(cls, prefixes):
        """Returns (size, number of named intervals, total number of intervals)
        for a set of areas of interest."""
        areas_of_interest = {}
        for (prefix_name, prefix) in prefixes.iteritems():
            prefix_label = prefix.get("Label", prefix_name)
            for (roi_name, roi_filename) in prefix.get("RegionsOfInterest", {}).iteritems():
                count, names, size = 0, set(), 0
                with open(roi_filename) as handle:
                    for line in handle:
                        bed = BEDRecord(line)
                        names.add(bed.name if len(bed) >= 4 else (bed.contig + "*"))
                        size += (bed.end - bed.start)
                        count += 1
                areas_of_interest[(prefix_name, roi_name)] = {"Size"       : size,
                                                              "NFeatures"  : len(names),
                                                              "NIntervals" : count,
                                                              "Genome"     : prefix["Name"],
                                                              "Name"       : roi_name,
                                                              "Label"      : "%s:%s" % (prefix_label, roi_name),
                                                              "Path"       : roi_filename}
        return areas_of_interest


    @classmethod
    def _stat_prefixes(cls, prefixes):
        """Returns (size, number of contigs) for a set of BWA prefix."""
        genomes = {}
        for prefix in prefixes:
            with open(prefixes[prefix]["Reference"] + ".fai") as table:
                lengths = [int(line.split()[1]) for line in table]

            labels = [prefix]
            if "Label" in prefixes[prefix]:
                labels.append(prefixes[prefix].get("Label"))

            for label in labels:
                if label not in genomes:
                    genomes[label] = {"Size" : 0, "NContigs" : 0}

                statistics = genomes[label]
                statistics["Size"] += sum(lengths)
                statistics["NContigs"] += len(lengths)

        if "mitochondrial" in genomes and "nuclear" in genomes:
            nucl = genomes["nuclear"]
            mito = genomes["mitochondrial"]

            genomes["endogenous"] = {"Size" : nucl["Size"] + mito["Size"],
                                     "NContigs" : nucl["NContigs"] + mito["NContigs"]}

        return genomes



def _measure_ordering(pair):
    measure = pair[0]
    key = measure.split("(")[0]
    return (__ORDERING[key], measure)


__ORDERING = {
    "lib_type"             :  00,
    "seq_reads_se"         :  10,
    "seq_trash_se"         :  20,
    "seq_trash_se_frac"    :  30,
    "seq_reads_pairs"      :  40,
    "seq_trash_pe_1"       :  50,
    "seq_trash_pe_1_frac"  :  60,
    "seq_trash_pe_2"       :  70,
    "seq_trash_pe_2_frac"  :  80,
    "seq_collapsed"        :  90,
    "seq_collapsed_frac"   : 100,
    "seq_retained_reads"   : 110,
    "seq_retained_nts"     : 120,
    "seq_retained_length"  : 130,


    "hits_raw"             : 140,
    "hits_raw_frac"        : 150,
    "hits_clonality"       : 160,
    "hits_unique"          : 170,
    "hits_unique_frac"     : 180,
    "hits_coverage"        : 190,
    "hits_length"          : 200,
    "ratio_reads"          : 210,
    "ratio_genome"         : 220,
    }
