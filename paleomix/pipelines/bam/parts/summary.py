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
import collections
import json
import math
import numbers
import os
import re
import sys

from paleomix.node import Node, NodeError
from paleomix.common.utilities import set_in, get_in
from paleomix.common.fileutils import move_file, reroot_path
from paleomix.tools.bam_stats.coverage import read_table as read_coverage_table
from paleomix.common.bedtools import BEDRecord
from paleomix.common.formats.fasta import FASTA

import paleomix.common.text as text


_PE_READS = frozenset(("Paired", "Singleton", "Collapsed", "CollapsedTruncated"))
_SE_READS = frozenset(("Single",))


class SummaryTableNode(Node):
    def __init__(
        self, config, makefile, target, cov_for_lanes, cov_for_libs, dependencies=()
    ):
        self._target = target.name
        self._output_file = os.path.join(config.destination, self._target + ".summary")
        self._prefixes = makefile["Prefixes"]

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
                        filenames = []
                        filetype = None

                        if lane.reads.stats:
                            filetype = "Raw"
                            filenames.append(lane.reads.stats)
                        elif lane.reads.validation:
                            filenames.extend(lane.reads.validation)

                            fileset = set(lane.reads.files)
                            if fileset & _PE_READS and fileset & _SE_READS:
                                filetype = "*"
                            elif fileset & _PE_READS:
                                filetype = "PE"
                            elif fileset & _SE_READS:
                                filetype = "SE"
                            else:
                                assert False, lane.reads.files
                                continue
                        else:
                            assert False

                        input_files.update(filenames)
                        self._in_raw_read[(sample.name, library.name, lane.name)] = (
                            filetype,
                            filenames,
                        )

        Node.__init__(
            self,
            description="writing summary to %s" % (self._output_file,),
            input_files=[_f for _f in input_files if _f],
            output_files=[self._output_file],
            dependencies=dependencies,
        )

    def _run(self, config, temp):
        rois = self._stat_areas_of_interest(self._prefixes)
        genomes = self._stat_prefixes(self._prefixes)
        with open(reroot_path(temp, self._output_file), "w") as table:
            table.write("# Command:\n")
            table.write("#     %s\n" % (" ".join(sys.argv)))
            table.write("#\n")
            table.write("# Directory:\n")
            table.write("#     %s\n" % (os.getcwd()))
            table.write("#\n")
            self._write_genomes(table, genomes)
            table.write("#\n")
            self._write_areas_of_interest(table, rois)
            table.write("#\n#\n")

            for roi in rois.values():
                genomes[roi["Label"]] = {"Size": roi["Size"]}
            self._write_tables(table, genomes)

    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._output_file), self._output_file)

    def _write_genomes(self, table, genomes):
        table.write("# Genomes:\n")
        rows = [["Name", "Contigs", "Size", "Prefix"]]
        for (_, prefix) in sorted(self._prefixes.items()):
            stats = genomes[prefix["Name"]]
            rows.append(
                (
                    prefix["Name"],
                    stats["NContigs"],
                    stats["Size"],
                    prefix["Path"],
                )
            )

        for line in text.padded_table(rows):
            table.write("#     %s\n" % (line,))

    def _write_areas_of_interest(self, table, rois):
        table.write("# Regions Of Interest:\n")
        rows = [["Genome", "ROI", "Size", "NFeatures", "NIntervals", "Path"]]
        for (_, roi) in sorted(rois.items()):
            rows.append(
                [
                    roi[key]
                    for key in (
                        "Genome",
                        "Name",
                        "Size",
                        "NFeatures",
                        "NIntervals",
                        "Path",
                    )
                ]
            )

        for line in text.padded_table(rows):
            table.write("#     %s\n" % (line,))

    def _write_tables(self, out, genomes):
        for row in self._build_tables(genomes):
            out.write("%s\n" % "\t".join(map(str, row)))

    def _build_tables(self, genomes):
        yield ["Target", "Sample", "Library", "Measure", "Value", "# Description"]

        for (target, samples) in sorted(
            self._read_tables(self._prefixes, genomes).items()
        ):
            for (sample, libraries) in sorted(samples.items()):
                for (library, prefixes) in sorted(libraries.items()):
                    ordered = (
                        [("reads", prefixes.pop("reads"))]
                        if "reads" in prefixes
                        else []
                    )
                    ordered.extend(sorted(prefixes.items()))

                    for (prefix, table) in ordered:
                        table.pop("hits_unique_nts(%s)" % prefix, None)

                        for (key, value) in sorted(
                            table.items(), key=_measure_ordering
                        ):
                            if isinstance(value, numbers.Number) and math.isnan(value):
                                value = "NA"

                            yield (
                                target,
                                sample,
                                library,
                                key,
                                value,
                                _get_comment(key),
                            )
                        yield ""
                    yield ""

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
                        merged_libraries[tblname] = self._merge_tables(
                            (merged_libraries.get(tblname, {}), subtable)
                        )
                        merged_samples[tblname] = self._merge_tables(
                            (merged_samples.get(tblname, {}), subtable)
                        )
                    libraries[library] = self._annotate_subtables(subtables, genomes)
                set_in(
                    table,
                    (target, sample, "*"),
                    self._annotate_subtables(merged_libraries, genomes),
                )
            set_in(
                table,
                (target, "*", "*"),
                self._annotate_subtables(merged_samples, genomes),
            )

        return table

    @classmethod
    def _annotate_subtables(cls, subtables, genomes):
        for (tblname, subtable) in subtables.items():
            if tblname == "reads":
                fractions = [
                    ("seq_trash_se", "seq_reads_se", "seq_trash_se_frac"),
                    ("seq_trash_pe_1", "seq_reads_pairs", "seq_trash_pe_1_frac"),
                    ("seq_trash_pe_2", "seq_reads_pairs", "seq_trash_pe_2_frac"),
                    ("seq_collapsed", "seq_reads_pairs", "seq_collapsed_frac"),
                    ("seq_retained_nts", "seq_retained_reads", "seq_retained_length"),
                ]

                for (numerator, denominator, measure) in fractions:
                    if (numerator in subtable) and (denominator in subtable):
                        value = subtable[numerator] / float(
                            subtable[denominator] or "NaN"
                        )
                        subtable[measure] = value
            else:
                total_hits = subtable["hits_raw(%s)" % tblname]
                total_nts = subtable["hits_unique_nts(%s)" % tblname]
                total_uniq = subtable["hits_unique(%s)" % tblname]
                total_reads = subtables.get("reads", {}).get("seq_retained_reads", 0)
                genome_size = genomes[tblname]["Size"]

                subtable["hits_raw_frac(%s)" % tblname] = total_hits / float(
                    total_reads or "NaN"
                )
                subtable["hits_unique_frac(%s)" % tblname] = total_uniq / float(
                    total_reads or "NaN"
                )
                subtable["hits_clonality(%s)" % tblname] = 1 - total_uniq / float(
                    total_hits or "NaN"
                )
                subtable["hits_length(%s)" % tblname] = total_nts / float(
                    total_uniq or "NaN"
                )
                subtable["hits_coverage(%s)" % tblname] = total_nts / float(
                    genome_size or "NaN"
                )

        return subtables

    def _read_reads_settings(self, table):
        for (
            (sample, library, barcode),
            (filetype, filenames),
        ) in self._in_raw_read.items():
            key = (self._target, sample, library, "reads", barcode)
            set_in(table, key, self._stat_read_settings(filetype, filenames))

        for (_, samples) in table.items():
            for (sample, libraries) in samples.items():
                for (library, prefixes) in libraries.items():
                    prefixes["reads"] = self._merge_tables(
                        list(prefixes["reads"].values())
                    )

        return table

    def _read_raw_bam_stats(self, table):
        for ((genome, target, sample, library), filenames) in self._in_raw_bams.items():
            key = (target, sample, library)
            hits, _ = self._read_coverage_tables(key, filenames)

            set_in(
                table, (target, sample, library, genome, "hits_raw(%s)" % genome), hits
            )

    def _read_lib_bam_stats(self, table):
        for ((genome, target, sample, library), filenames) in self._in_lib_bams.items():
            key = (target, sample, library)
            hits, nts = self._read_coverage_tables(key, filenames)

            set_in(
                table,
                (target, sample, library, genome, "hits_unique(%s)" % genome),
                hits,
            )
            set_in(
                table,
                (target, sample, library, genome, "hits_unique_nts(%s)" % genome),
                nts,
            )

    @classmethod
    def _read_coverage_tables(cls, key, filenames):
        hits = nts = 0
        for filename in filenames:
            subtable = {}
            read_coverage_table(subtable, filename)
            contigtables = get_in(subtable, key)

            if contigtables is None:
                raise NodeError(
                    "Error reading table %r; row not found:"
                    "\n   %s\n\nIf files have been renamed "
                    "during the run, then please remove this file "
                    "in that it may be re-generated.\nHowever, "
                    "note that read-group tags in the BAM files "
                    "may not be correct!" % (filename, "   ".join(key))
                )

            for contigtable in contigtables.values():
                hits += contigtable["Hits"]
                nts += contigtable["M"]
        return hits, nts

    @classmethod
    def _merge_tables(cls, tables):
        merged = {}
        for table in tables:
            for (measure, value) in table.items():
                if not isinstance(value, numbers.Number):
                    other = merged.get(measure, value)
                    merged[measure] = value if (value == other) else "*"
                else:
                    merged[measure] = value + merged.get(measure, 0)
        return merged

    @classmethod
    def _stat_read_settings(cls, filetype, filenames):
        assert filetype in ("Raw", "SE", "PE", "*", "BAM"), filetype

        if filetype != "Raw":
            nan = float("nan")
            stats = {
                "lib_type": filetype,
                "seq_retained_nts": 0,
                "seq_retained_reads": 0,
            }

            if filetype in ("SE", "*"):
                stats["seq_reads_se"] = nan
                stats["seq_trash_se"] = nan

            if filetype in ("PE", "*"):
                stats["seq_reads_pairs"] = nan
                stats["seq_trash_pe_1"] = nan
                stats["seq_trash_pe_2"] = nan
                stats["seq_collapsed"] = 0

            for filename in filenames:
                with open(filename) as handle:
                    for key, value in json.load(handle).items():
                        if key in stats:
                            stats[key] += int(value)

            return stats

        assert len(filenames) == 1, filenames
        with open(filenames[0]) as settings_file:
            settings = settings_file.read()

            def _re_search(regexp, default=None):
                match = re.search(regexp, settings)
                if not match:
                    if default is not None:
                        return default
                    raise KeyError(
                        "Could not find match with RegExp %s in file '%s'"
                        % (repr(regexp), filename)
                    )

                return int(match.groups()[0])

            if "Paired end mode" in settings or "paired-end reads" in settings:
                return {
                    "lib_type": "PE",
                    "seq_reads_pairs": _re_search("number of read pairs: ([0-9]+)"),
                    "seq_trash_pe_1": _re_search("discarded mate 1 reads: ([0-9]+)"),
                    "seq_trash_pe_2": _re_search("discarded mate 2 reads: ([0-9]+)"),
                    "seq_retained_nts": _re_search("retained nucleotides: ([0-9]+)"),
                    "seq_retained_reads": _re_search("retained reads: ([0-9]+)"),
                    "seq_collapsed": _re_search(
                        "of (?:full-length )?collapsed pairs: ([0-9]+)", 0
                    )
                    + _re_search("of truncated collapsed pairs: ([0-9]+)", 0),
                }
            elif "Single end mode" in settings or "single-end reads" in settings:
                return {
                    "lib_type": "SE",
                    "seq_reads_se": _re_search(
                        "number of (?:reads|read pairs): ([0-9]+)"
                    ),
                    "seq_trash_se": _re_search("discarded mate 1 reads: ([0-9]+)"),
                    "seq_retained_nts": _re_search("retained nucleotides: ([0-9]+)"),
                    "seq_retained_reads": _re_search("retained reads: ([0-9]+)"),
                }
            else:
                assert False, filename

    @classmethod
    def _stat_areas_of_interest(cls, prefixes):
        """Returns (size, number of named intervals, total number of intervals)
        for a set of areas of interest."""
        areas_of_interest = {}
        for (prefix_name, prefix) in prefixes.items():
            for (roi_name, roi_filename) in prefix.get("RegionsOfInterest", {}).items():
                count, names, size = 0, set(), 0
                with open(roi_filename) as handle:
                    for line in handle:
                        bed = BEDRecord(line)
                        names.add(bed.name if len(bed) >= 4 else (bed.contig + "*"))
                        size += bed.end - bed.start
                        count += 1
                areas_of_interest[(prefix_name, roi_name)] = {
                    "Size": size,
                    "NFeatures": len(names),
                    "NIntervals": count,
                    "Genome": prefix["Name"],
                    "Name": roi_name,
                    "Path": roi_filename,
                }
        return areas_of_interest

    @classmethod
    def _stat_prefixes(cls, prefixes):
        """Returns (size, number of contigs) for a set of BWA prefix."""
        genomes = {}
        for prefix in prefixes:
            contigs = FASTA.index_and_collect_contigs(prefixes[prefix]["Reference"])

            genomes[prefix] = {
                "Size": sum(contigs.values()),
                "NContigs": len(contigs),
            }

        return genomes


def _measure_ordering(pair):
    measure = pair[0]
    key = measure.split("(")[0]
    return (_ORDERING[key], measure)


def _get_comment(key):
    comment = _COMMENTS.get(key)
    if comment is None:
        comment = _COMMENTS[key.split("(", 1)[0]]

    return comment


_ORDERING = {
    "lib_type": 00,
    "seq_reads_se": 10,
    "seq_trash_se": 20,
    "seq_trash_se_frac": 30,
    "seq_reads_pairs": 40,
    "seq_trash_pe_1": 50,
    "seq_trash_pe_1_frac": 60,
    "seq_trash_pe_2": 70,
    "seq_trash_pe_2_frac": 80,
    "seq_collapsed": 90,
    "seq_collapsed_frac": 100,
    "seq_retained_reads": 110,
    "seq_retained_nts": 120,
    "seq_retained_length": 130,
    "hits_raw": 140,
    "hits_raw_frac": 150,
    "hits_clonality": 160,
    "hits_unique": 170,
    "hits_unique_frac": 180,
    "hits_coverage": 190,
    "hits_length": 200,
    "ratio_reads": 210,
    "ratio_genome": 220,
}

_COMMENTS = {
    "hits_clonality": "# Fraction of hits that were PCR duplicates",
    "hits_coverage": "# Estimated coverage from unique hits",
    "hits_length": "# Average number of aligned bases per unique hit",
    "hits_raw_frac": "# Total number of hits vs. total number of reads retained",
    "hits_raw": "# Total number of hits (prior to PCR duplicate filtering)",
    "hits_unique_frac": "# Total number of unique hits vs. total number of reads retained",
    "hits_unique": "# Total number of hits (excluding any PCR duplicates)",
    "lib_type": "# SE, PE, or * (for both)",
    "ratio_genome(mito,nuc)": "# Ratio of NTs of unique hits corrected by genome sizes: (NTs(mito) / NTs(nuc)) / (Size(mito) / (2 * Size(nuc)))",
    "ratio_genome(nuc,mito)": "# Ratio of NTs of unique hits corrected by genome sizes: (NTs(nuc) / NTs(mito)) / ((2 * Size(nuc)) / Size(mito))",
    "ratio_reads(nuc,mito)": "# Ratio of unique hits: Hits(nuc) / H(mito)",
    "seq_collapsed_frac": "# Fraction of PE pairs collapsed into one read",
    "seq_collapsed": "# Total number of pairs collapsed into one read",
    "seq_reads_pairs": "# Total number of pairs",
    "seq_reads_se": "# Total number of single-ended reads",
    "seq_retained_length": "# Average number of NTs in retained reads",
    "seq_retained_nts": "# Total number of NTs in retained reads",
    "seq_retained_reads": "# Total number of retained reads",
    "seq_trash_pe_1_frac": "# Fraction of PE mate 1 reads trashed",
    "seq_trash_pe_1": "# Total number of reads",
    "seq_trash_pe_2_frac": "# Fraction of PE mate 2 reads trashed",
    "seq_trash_pe_2": "# Total number of reads",
    "seq_trash_se_frac": "# Fraction of SE reads trashed",
    "seq_trash_se": "# Total number of trashed reads",
}
