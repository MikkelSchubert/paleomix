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
import subprocess
import collections

import pysam

from pypeline.node import Node
from pypeline.common.utilities import safe_coerce_to_tuple, set_in, get_in
from pypeline.common.fileutils import swap_ext, move_file, reroot_path
from pypeline.nodes.coverage import read_table as read_coverage_table

import pypeline.common.text as text
import pypeline.tools.bam_pipeline as common




class SummaryTableNode(Node):
    def __init__(self, config, makefile, prefixes, target, samples, records, dependencies = ()):
        self._target        = target
        self._output_file   = os.path.join(config.destination, self._target + ".summary")
        self._prefixes      = prefixes
        self._records       = {target : samples}
        self._makefile      = makefile

        self._in_raw_read = collections.defaultdict(list)
        for (sample, libraries) in samples.iteritems():
            for (library, barcodes) in libraries.iteritems():
                for (barcode, record) in barcodes.iteritems():
                    self._in_raw_read[(sample, library, barcode)] = None
                    if "Raw" in record["Reads"]:
                        self._in_raw_read[(sample, library, barcode)] = \
                            os.path.join(config.destination, target, "reads", sample, library, barcode, "reads.settings")

        self._in_raw_bams = collections.defaultdict(list)
        self._in_lib_bams = collections.defaultdict(list)
        for genome in prefixes:
            label = prefixes[genome].get("Label") or genome
            for record in records[genome]:
                for (key, dd) in (("LaneCoverage", self._in_raw_bams), ("LibCoverage", self._in_lib_bams)):
                    for node in record[key]:
                        dd[(label, target, record["Sample"], record["Library"])].extend(node.output_files)

        
        input_files = self._in_raw_read.values() \
            + sum(self._in_raw_bams.values(), []) \
            + sum(self._in_lib_bams.values(), [])

        Node.__init__(self,
                      description  = "<Summary: %s>" % self._output_file,
                      input_files  = filter(None, input_files),
                      output_files = [self._output_file],
                      dependencies = dependencies)


    def _run(self, config, temp):
        genomes = self._stat_bwa_prefixes(self._prefixes)
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
            table.write("#\n#\n")
            self._write_tables(table, genomes)


    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._output_file), self._output_file)


    def _write_genomes(self, table, genomes):
        table.write("# Genomes:\n")
        rows = [["Name", "Label", "Contigs", "Size", "Prefix"]]
        for (_, prefix) in sorted(self._prefixes.items()):
            stats = genomes[prefix.get("Label") or prefix["Name"]]
            rows.append((prefix["Name"], (prefix.get("Label") or "-"), stats["NContigs"], stats["Size"], prefix["Path"]))

        for line in text.padded_table(rows):
            table.write("#     %s\n" % (line,))


    def _write_tables(self, out, genomes):
        rows = [["Target", "Sample", "Library", "Measure", "Value", ""]]
        for (target, samples) in sorted(self._read_tables(self._prefixes, genomes).iteritems()):
            for (sample, libraries) in sorted(samples.iteritems()):
                for (library, prefixes) in sorted(libraries.iteritems()):
                    ordered = [("reads", prefixes.pop("reads"))]
                    ordered.extend(sorted(prefixes.items()))
                    
                    for (prefix, table) in ordered:
                        table.pop("hits_unique_nts(%s)" % prefix, None)

                        for (key, (value, comment)) in sorted(table.iteritems(), key = _measure_ordering):
                            if isinstance(value, numbers.Number) and math.isnan(value):
                                value = "NA"
                            rows.append((target, sample, library, key, value, comment))
                        rows.append([])
                    rows.append([])
 
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
                total_reads = subtables["reads"]["seq_retained_reads"][0]
                                
                subtable["hits_raw_frac(%s)" % tblname] = (total_hits / float(total_reads), "# Total number of hits vs. total number of reads")
                subtable["hits_unique_frac(%s)" % tblname] = (total_uniq / float(total_reads), "# Total number of unique hits vs. total number of reads")
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
            subtable = {}
            for filename in filenames:
                read_coverage_table(subtable, filename)
            key = (target, sample, library)

            hits = 0
            for contigtable in get_in(subtable, key).itervalues():
                hits += contigtable["Hits"]
            
            value = (hits, "# Total number of hits (prior to PCR duplicate filtering)")
            set_in(table, (target, sample, library, genome, "hits_raw(%s)" % genome), value)


    def _read_lib_bam_stats(self, table):
        for ((genome, target, sample, library), filenames) in self._in_lib_bams.iteritems():
            subtable = {}
            for filename in filenames:
                read_coverage_table(subtable, filename)
            key = (target, sample, library)

            hits = nts = 0
            for contigtable in get_in(subtable, key).itervalues():
                hits += contigtable["Hits"]
                nts += contigtable["M"]
            
            value = (hits, "# Total number of hits (excluding any PCR duplicates)")
            set_in(table, (target, sample, library, genome, "hits_unique(%s)" % genome), value)
            set_in(table, (target, sample, library, genome, "hits_unique_nts(%s)" % genome), (nts, None))


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
        if filename is None:
            # FIXME: A better solution is required when adding support pre-trimmed reads ...
            return {
                "lib_type"           : ("*", "# SE, PE, or * (for both)"),
                "seq_reads_se"       : (float("nan"),  "# Total number of single-ended reads"),
                "seq_trash_se"       : (float("nan"),  "# Total number of trashed reads"),
                "seq_retained_nts"   : (float("nan"),  "# Total number of NTs in retained reads"),
                "seq_retained_reads" : (float("nan"),  "# Total number of retained reads")
                }
        
        with open(filename) as settings_file:
            settings = settings_file.read()
            if "Paired end mode" in settings:
                return {
                    "lib_type"            : ("PE", "# SE, PE, or * (for both)"),
                    "seq_reads_pairs"     : (int(re.search("read pairs: ([0-9]+)",             settings).groups()[0]),  "# Total number of reads"),
                    "seq_trash_pe_1"      : (int(re.search("discarded mate 1 reads: ([0-9]+)", settings).groups()[0]),  "# Total number of reads"),
                    "seq_trash_pe_2"      : (int(re.search("discarded mate 2 reads: ([0-9]+)", settings).groups()[0]),  "# Total number of reads"),
                    "seq_retained_nts"    : (int(re.search("retained nucleotides: ([0-9]+)",   settings).groups()[0]),  "# Total number of NTs in retained reads"),
                    "seq_retained_reads"  : (int(re.search("retained reads: ([0-9]+)",         settings).groups()[0]),  "# Total number of retained reads"),
                    "seq_collapsed"       : (int(re.search("well aligned pairs: ([0-9]+)",     settings).groups()[0]),  "# Total number of pairs collapsed into one read"),
                    }
            elif "Single end mode" in settings:
                return {
                    "lib_type"            : ("SE", "# SE, PE, or * (for both)"),
                    "seq_reads_se"        : (int(re.search("read pairs: ([0-9]+)",             settings).groups()[0]),  "# Total number of single-ended reads"),
                    "seq_trash_se"        : (int(re.search("discarded mate 1 reads: ([0-9]+)", settings).groups()[0]),  "# Total number of trashed reads"),
                    "seq_retained_nts"    : (int(re.search("retained nucleotides: ([0-9]+)",   settings).groups()[0]),  "# Total number of NTs in retained reads"),
                    "seq_retained_reads"  : (int(re.search("retained reads: ([0-9]+)",         settings).groups()[0]),  "# Total number of retained reads"),
                    }
            else:
                assert False


    @classmethod
    def _stat_bwa_prefixes(cls, prefixes):
        """Returns (size, number of contigs) for a set of BWA prefix."""
        genomes = {}
        for prefix in prefixes:
            label = prefixes[prefix].get("Label") or prefix
            with open(prefixes[prefix]["Path"] + ".ann") as table:
                genomes[label] = dict(zip(("Size", "NContigs"), map(int, table.readline().strip().split())[:2]))
        
        if "mitochondrial" in genomes and "nuclear" in genomes:
            nucl = genomes["nuclear"]
            mito = genomes["mitochondrial"]

            genomes["endogenous"] = {"Size" : nucl["Size"] + mito["Size"],
                                     "NContigs" : nucl["NContigs"] + mito["NContigs"]}

        return genomes



def _measure_ordering((measure, _)):
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
