import os
import re
import sys
import numbers
import subprocess
import collections

import pysam

from pypeline.common.fileutils import swap_ext, move_file, reroot_path

import pypeline.common.text as text
import pypeline.tools.alignment_common as common
from pypeline.node import Node



class SummaryTableNode(Node):
    def __init__(self, config, records, dependencies = ()):
        self._records = records
        self._target = set(record["Name"] for record in records)
        assert len(self._target) == 1, self._target
        self._target = self._target.pop()
        self._output_file = os.path.join("results", self._target + ".summary")

        input_files = set()
        for prefix in config.bwa_prefix:
            for record in records:
                library_prefix = os.path.join(common.paths.sample_path(record, prefix), record["Library"])

                input_files.add(library_prefix + ".unaligned.bam")
                input_files.add(library_prefix + ".unaligned.markdup.bam")
                if record["files"]["R2"]:
                    input_files.add(library_prefix + ".aligned.bam")
                    input_files.add(library_prefix + ".aligned.kirdup.bam")
        input_files |= set(swap_ext(filename, ".bai") for filename in input_files)
        

        Node.__init__(self,
                      description  = "<Summary: %s>" % self._output_file,
                      input_files  = input_files,
                      output_files = [self._output_file],
                      dependencies = dependencies)
        


    def _run(self, config, temp):
        with open(reroot_path(temp, self._output_file), "w") as table:
            table.write("# Command:\n")
            table.write("#     %s\n" % (" ".join(sys.argv)),)
            table.write("#\n")
            table.write("# Directory:\n")
            table.write("#     %s\n" % (os.getcwd()),)
            table.write("#\n")
            self._write_records(table, self._records)
            table.write("#\n")
            self._write_genomes(table, config)
            table.write("#\n#\n")
            self._write_table(table, config, self._records)


    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._output_file), self._output_file)


    def _write_records(self, table, records):
        rows = [["Name", "Sample", "Library", "Barcode", "Platform", "Path"]]
        for record in records:
            rows.append([record[field] for field in rows[0]])

        table.write("# Records:\n")
        for line in text.padded_table(rows):
            table.write("#     %s\n" % (line,))        


    def _write_genomes(self, table, config):
        table.write("# Genomes:\n")
        rows = [["Name", "Label", "Contigs", "Size", "Prefix"]]
        for prefix in config.bwa_prefix:
            label = _prefix_to_label(config, prefix, only_endogenous = True)
            size, contigs = self._stat_bwa_prefix(prefix)
            rows.append((common.paths.prefix_to_filename(prefix), label, contigs, size, prefix))

        for line in text.padded_table(rows):
            table.write("#     %s\n" % (line,))

    
    def _write_table(self, table, config, records):
        annotated_table = self._annotate_summary_table(config, self._build_summary_table(config, records))

        table_rows = [["Name", "Sample", "Library", "Measure", "Value", ""]]
        for ((sample, library), subtable) in sorted(annotated_table.items()):
            for (measure, (value, comment)) in sorted(subtable["seq"].iteritems(), key = _measure_ordering):
                table_rows.append([self._target, sample, library, measure, value, comment])
            table_rows.extend("##")

            for (genome, hit_table) in sorted(subtable["hits"].items()):
                for (measure, (value, comment)) in sorted(hit_table.iteritems(), key = _measure_ordering):
                    table_rows.append([self._target, sample, library, measure, value, comment])
                table_rows.append("#")
            table_rows.extend("##")

        while table_rows[-1] == "#":
            table_rows.pop()               
        
        for line in text.padded_table(table_rows):
            table.write("%s\n" % (line,))


    @classmethod
    def _annotate_summary_table(cls, config, table):
        fractions = [("seq_trash_se",   "seq_reads_se",    "seq_trash_se_frac",   "# Fraction of SE reads trashed"),
                     ("seq_trash_pe_1", "seq_reads_pairs", "seq_trash_pe_1_frac", "# Fraction of PE mate 1 reads trashed"), 
                     ("seq_trash_pe_2", "seq_reads_pairs", "seq_trash_pe_2_frac", "# Fraction of PE mate 2 reads trashed"), 
                     ("seq_collapsed",  "seq_reads_pairs", "seq_collapsed_frac",  "# Fraction of PE pairs collapsed into one read")]

        genomes = {}
        for prefix in config.bwa_prefix:
            label = _prefix_to_label(config, prefix, only_endogenous = False)
            genomes[label] = cls._stat_bwa_prefix(prefix)[0]


        for ((sample, library), level) in table.iteritems():
            sequences = level["seq"]
            for (numerator, denominator, measure, comment) in fractions:
                if (numerator in sequences) and (denominator in sequences):
                    value = float(sequences[numerator][0]) / sequences[denominator][0]
                    sequences[measure] = (value, comment)
                
            total_reads = 0
            total_reads += sequences.get("seq_reads_se",    (0,))[0] - sequences.get("seq_trash_se", (0,))[0]
            total_reads += sequences.get("seq_reads_pairs", (0,))[0] - sequences.get("seq_trash_pe_1", (0,))[0]
            total_reads += sequences.get("seq_reads_pairs", (0,))[0] - sequences.get("seq_trash_pe_2", (0,))[0]
            sequences["seq_reads_total"] = (total_reads, )

            average_nts = sequences.get("seq_nt_total", (0,))[0] / float(total_reads)
            sequences["seq_nt_average"] = (average_nts, "# Average number of NTs in retained reads")

            
            hits = level["hits"]
            for (genome, hit_table) in hits.iteritems():
                total_hits            = hit_table["hits_raw(%s)" % genome][0]
                total_hits_unique     = hit_table["hits_unique(%s)" % genome][0]
                total_hits_unique_nts = hit_table["hits_unique_nts(%s)" % genome][0]

                hit_table["hits_raw_frac(%s)" % genome]       = (float(total_hits) / total_reads,                  "# Fraction of reads succesfully mapped")
                hit_table["hits_coverage(%s)" % genome]       = (float(total_hits_unique_nts) / genomes[genome],   "# Estimated coverage from unique hits")
                if total_hits:
                    hit_table["hits_clonality(%s)"  % genome] = (1 - float(total_hits_unique) / total_hits,        "# Fraction of hits that were PCR duplicates")
                    hit_table["hits_length(%s)" % genome]     = (float(total_hits_unique_nts) / total_hits_unique, "# Average length of unique hits")
                else:
                    hit_table["hits_clonality(%s)"  % genome] = ("NA", "# Fraction of hits that were PCR duplicates")
                    hit_table["hits_length(%s)" % genome]     = ("NA", "# Average length of unique hits")



            if ("mito" in level["hits"]) and ("nuclear" in level["hits"]):
                mito_size, _ = cls._stat_bwa_prefix(config.bwa_prefix_mito)
                nuc_size,  _ = cls._stat_bwa_prefix(config.bwa_prefix_nuclear)
                genome_size  = mito_size + nuc_size

                total_hits            = hits["mito"]["hits_raw(mito)"][0]         + hits["nuclear"]["hits_raw(nuclear)"][0]
                total_hits_unique     = hits["mito"]["hits_unique(mito)"][0]      + hits["nuclear"]["hits_unique(nuclear)"][0]
                total_hits_unique_nts = hits["mito"]["hits_unique_nts(mito)"][0]  + hits["nuclear"]["hits_unique_nts(nuclear)"][0]
                total_reads           = sequences["seq_reads_total"][0]

                clonality = "NA" if not total_hits else (1 - float(total_hits_unique) / total_hits)
                coverage  = float(total_hits_unique_nts) / (genomes["mito"] + genomes["nuclear"])
                length    = float(total_hits_unique_nts) / total_hits_unique

                ratio, genome_rt, genome_iv = "NA", "NA", "NA"
                if hits["mito"]["hits_unique(mito)"][0]:
                    ratio     = float(hits["nuclear"]["hits_unique(nuclear)"][0]) / hits["mito"]["hits_unique(mito)"][0]
                    genome_rt = ratio * float(genomes["mito"]) / genomes["nuclear"]
                    genome_iv = genome_rt ** -1

                hits["endogenous"] = {
                    "hits_raw(endogenous)"       : (total_hits,                       "# Total number of hits against the nuclear and mitochondrial genome"),
                    "hits_raw_frac(endogenous)"  : (float(total_hits) / total_reads,  "# Fraction of reads that mapped against the nuclear and mitochondrial genome"),
                    "hits_unique(endogenous)"    : (total_hits_unique,                "# Total number of unique reads (PRC duplicates removed)"),
                    "hits_clonality(endogenous)" : (clonality,                        "# Fraction of hits that were PCR duplicates"),
                    "hits_coverage(endogenous)"  : (coverage,                         "# Estimated coverage from unique hits"),
                    "hits_length(endogenous)"    : (length,                           "# Average length of unique hits"),
                    "ratio_reads(nuc,mito)"      : (ratio,                            "# Ratio of unique hits"),
                    "ratio_genome(nuc,mito)"     : (genome_rt,                        "# Ratio of unique hits (corrected by genome sizes)"),
                    "ratio_genome(mito,nuc)"     : (genome_iv,                        "# Ratio of unique hits (corrected by genome sizes)"),
                    }


            # Cleanup 
            del sequences["seq_nt_total"]
            del sequences["seq_reads_total"]
            for prefix in config.bwa_prefix:
                label = _prefix_to_label(config, prefix, only_endogenous = False)
                del hits[label]["hits_unique_nts(%s)" % label]
                

        return table
                

    @classmethod
    def _build_summary_table(cls, config, records):
        sequences = collections.defaultdict(list)
        for record in records:
            key = (record["Sample"], record["Library"])
            sequences[key].append(cls._stat_read_settings(record))
            
        hits = {}
        for prefix in config.bwa_prefix:
            label = _prefix_to_label(config, prefix, only_endogenous = False)
            dd = hits[label] = collections.defaultdict(list)
            
            for record in records:
                key = (record["Sample"], record["Library"])
                dd[key] = [cls._stat_hits(config, prefix, record)]

        keyfuncs = (lambda (sample, library): (sample, library),
                    lambda (sample, _): (sample, "*"),
                    lambda _: ("*", "*"))
        
        levels = collections.defaultdict(lambda: {"hits" : {}})
        for keyfunc in keyfuncs:
            for (key, value) in cls._merge_tables(sequences, keyfunc).iteritems():
                levels[key]["seq"] = value

            for prefix in config.bwa_prefix:
                label = _prefix_to_label(config, prefix, only_endogenous = False)
                for (key, value) in cls._merge_tables(hits[label], keyfunc).iteritems():
                    levels[key]["hits"][label] = value

        return levels


    @classmethod
    def _stat_hits(cls, config, prefix, record):
        label = _prefix_to_label(config, prefix, only_endogenous = False)
        library_prefix = common.paths.library_path(record, prefix)
        
        files_raw    = [library_prefix + ".unaligned.bam"]
        files_unique = [library_prefix + ".unaligned.markdup.bam"]
        if record["files"]["R2"]:
            files_raw    += [library_prefix + ".aligned.bam"]
            files_unique += [library_prefix + ".aligned.kirdup.bam"]

        raw_hits    = cls._stat_bam_index(files_raw)
        unique_hits = cls._stat_bam_index(files_unique)
        unique_nts  = cls._stat_bam_files(files_unique)
 
        return {
            "hits_raw(%s)"        % label : (raw_hits,    "# Total number of hits (including PCR duplicates)"),
            "hits_unique(%s)"     % label : (unique_hits, "# Total number of unique hits (excluding PCR duplicates)"),
            "hits_unique_nts(%s)" % label : (unique_nts,  ""),
            }


    @classmethod
    def _stat_bam_index(cls, filenames):
        hits = 0
        for filename in filenames:
            command = ["samtools", "idxstats", filename]
            proc = subprocess.Popen(command, stdout = subprocess.PIPE)
            for line in proc.stdout:
                fields = line.split("\t")
                if fields[0] != "*":
                    hits += int(fields[2])
                    
            if proc.wait():
                raise RuntimeError("Error running samtools idxstats on '%s': %i" % (filename, proc.returncode))

        return hits
            


    @classmethod
    def _stat_bam_files(cls, filenames):
        nucleotides = 0
        for filename in filenames:
            bamfile = pysam.Samfile(filename)
            try:
                nucleotides += sum(read.qlen for read in bamfile)
            finally:
                bamfile.close()
        
        return nucleotides
            


    @classmethod
    def _stat_bwa_prefix(cls, prefix):
        """Returns (size, number of contigs) for a given BWA prefix."""        
        with open(prefix + ".ann") as table:
            return map(int, table.readline().strip().split())[:2]

    
    @classmethod
    def _merge_tables(cls, stats, keyfunc):
        merged_tables = collections.defaultdict(dict)
        for (key, tables) in stats.iteritems():
            sample, library = keyfunc(key)

            merged = merged_tables[(sample, library)]
            for table in tables:
                for (measure, (value, comment)) in table.iteritems():
                    if not isinstance(value, numbers.Number):
                        other, _ = merged.get(measure, (value, None))
                        merged[measure] = (value if (value == other) else "*", comment)
                    else:
                        other, _ = merged.get(measure, (0, None))
                        merged[measure] = (value + other, comment)

        return merged_tables
        
 
        
    @classmethod
    def _stat_read_settings(cls, record):
        # FIXME: Add as input!
        with open(os.path.join(common.paths.full_path(record), "reads.settings")) as settings:
            text = settings.read()
            if "Paired end mode" in text:
                return {
                    "lib_type"        : ("PE", "# SE, PE, or * (for both)"),
                    "seq_reads_pairs"   : (int(re.search("read pairs: ([0-9]+)", text).groups()[0]),  "# Total number of reads"),
                    "seq_trash_pe_1"    : (int(re.search("discarded mate 1 reads: ([0-9]+)", text).groups()[0]),  "# Total number of reads"),
                    "seq_trash_pe_2"    : (int(re.search("discarded mate 2 reads: ([0-9]+)", text).groups()[0]),  "# Total number of reads"),
                    "seq_nt_total"      : (int(re.search("retained nucleotides: ([0-9]+)", text).groups()[0]),  "# Total number of NTs in retained reads"),
                    "seq_collapsed"     : (int(re.search("well aligned pairs: ([0-9]+)", text).groups()[0]),  "# Total number of pairs collapsed into one read"),
                    }
            elif "Single end mode" in text:
                return {
                    "lib_type"          : ("SE", "# SE, PE, or * (for both)"),
                    "seq_reads_se"      : (int(re.search("read pairs: ([0-9]+)", text).groups()[0]),  "# Total number of single-ended reads"),
                    "seq_trash_se"      : (int(re.search("discarded mate 1 reads: ([0-9]+)", text).groups()[0]),  "# Total number of reads"),
                    "seq_nt_total"      : (int(re.search("retained nucleotides: ([0-9]+)", text).groups()[0]),  "# Total number of NTs in retained reads")
                    }

            else:
                assert False




def _prefix_to_label(config, prefix, only_endogenous):
    if (prefix == config.bwa_prefix_mito):
        return "mito"
    elif (prefix == config.bwa_prefix_nuclear):
        return "nuclear"
    elif only_endogenous:
        return "-"

    return common.paths.prefix_to_filename(prefix)

    


def _measure_ordering((measure, _)):
    key = measure.split("(")[0]
    return (__ORDERING[key], measure)


__ORDERING = {
    "lib_type"             :  0,
    "seq_reads_se"         :  1,
    "seq_trash_se"         :  2,
    "seq_trash_se_frac"    :  3,
    "seq_reads_pairs"      :  4,
    "seq_trash_pe_1"       :  5,
    "seq_trash_pe_1_frac"  :  6,
    "seq_trash_pe_2"       :  7,
    "seq_trash_pe_2_frac"  :  8,
    "seq_collapsed"        :  9,
    "seq_collapsed_frac"   : 10,
    "seq_nt_average"       : 11,

    "hits_raw"             : 12,
    "hits_raw_frac"        : 13,
    "hits_clonality"       : 14,
    "hits_unique"          : 15,
    "hits_coverage"        : 16,
    "hits_length"          : 17,
    "ratio_reads"          : 18,
    "ratio_genome"         : 19,
    }
