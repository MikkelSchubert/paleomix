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
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
# Converts a GTF file to a set of BED6 files, one for each
# feature in the GTF file (CDS, Exon, ...). Also generates a list
# of introns, and UTRs based on sequences in the GTF.
#
from __future__ import with_statement
from __future__ import print_function

import sys
from argparse import ArgumentParser

import pysam

from paleomix.common.fileutils import open_ro
from paleomix.common.utilities import set_in, get_in

import paleomix.common.text as text


###############################################################################
###############################################################################
# Functions used for GTF parsing and filtering

def filter_gtf_record(gtf):
    return gtf.feature not in ("exon", "CDS")


def update_gtf_table(table, gtf, scaffolds, contig_prefix):
    # Workaround for bug in Pysam, which mis-parses individual properties
    # (e.g. exon_number) if these are not quoted. This does not apply to
    # asDict, which uses a different parsing implementation (v0.7.8).
    properties = gtf.asDict()

    gene_type = properties.get("gene_biotype")
    if gene_type is None:
        gene_type = properties.get("gene_type", "unknown_genetype")

    keys = (gene_type,
            properties["gene_id"],
            properties["transcript_id"],
            int(properties["exon_number"]),
            gtf.feature)

    record = {"contig": contig_prefix + gtf.contig,
              "start": gtf.start,
              # In pysam, 'end' equals the past-the-end position
              "end": gtf.end - 1,
              "strand": gtf.strand,
              "feature": gtf.feature,
              "transcript": properties["transcript_id"]}

    if record["contig"] in scaffolds:
        contig = scaffolds[record["contig"]]
        record["contig"] = contig["chrom"]
        record["start"] += int(contig["chromStart"])
        record["end"] += int(contig["chromStart"])

    assert not get_in(table, keys), keys
    set_in(table, keys, record)


def read_gtf(lines, scaffolds, contig_prefix):
    table = {}  # gene_id -> transcript_id -> exon_number -> feature -> [items]
    for gtf in text.parse_lines(lines, pysam.asGTF()):
        if not filter_gtf_record(gtf):
            update_gtf_table(table, gtf, scaffolds, contig_prefix)

    return table


def read_scaffolds(filename):
    scaffolds = {}
    with open(filename) as handle:
        header = handle.readline().strip("#\n\r").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\r\n").split("\t")))
            scaffolds[row["contig"]] = row
    return scaffolds


###############################################################################
###############################################################################

def get_introns(exons):
    lst = [exons[n]["exon"] for n in sorted(exons)]
    if lst[0]["strand"] == "-":
        # Exon numbers follow the read direction, and the coordinates are
        # therefore descending for regions on the negative strand. Below we
        # assume that the coordinates are ascending, so reorder the list.
        lst = lst[::-1]

    introns = []
    for (record_a, record_b) in zip(lst, lst[1:]):
        if record_a["end"] == record_b["start"] - 1:
            # Intron has been lost?
            continue

        record = dict(record_a)
        record.update(feature="intron",
                      start=record_a["end"] + 1,
                      end=record_b["start"] - 1)
        assert record["start"] <= record["end"], lst

        introns.append(record)

    return introns


def split_exon(exon, cds):
    """Takes an exon and a CDS, and returns a map of regions for each
       feature (UTR5/3, CDS) that may be inferred from the arguments.
       Note that the CDS is simply returned as is, to simplify
       downstream handling of these features."""
    results = [cds]

    if exon["start"] < cds["start"]:
        utr = dict(exon)
        utr.update(end=cds["start"] - 1,
                   feature=(exon["strand"] == "+" and "UTR5" or "UTR3"))
        results.append(utr)

    if exon["end"] > cds["end"]:
        utr = dict(exon)
        utr.update(start=cds["end"] + 1,
                   feature=(exon["strand"] == "+" and "UTR3" or "UTR5"))
        results.append(utr)

    return results


def split_exons(exons, func):
    # By looping over the list sorted by exon-number, we can easily
    # determine whether or not we are dealing with a 5' or 3' UTR.
    seen_cds = False
    for (_, exon) in sorted(exons.iteritems()):
        if "CDS" in exon:
            seen_cds = True
            cds, exon = exon["CDS"], exon["exon"]

            func(split_exon(exon, cds))
        else:
            utr = dict(exon["exon"])
            utr.update(feature=(seen_cds and "UTR3" or "UTR5"))

            func([utr])


def select_transcripts(options, transcripts, protein_coding):
    """Returns the largest transcript, preferring well formed
    transcripts (len(CDS) % 3 == 0) if the gene is protein coding."""
    if options.keep_all_transcripts and options.keep_malformed_proteins:
        return transcripts.itervalues()

    selection = []
    for transcript in transcripts.itervalues():
        well_formed = True
        exon_len = cds_len = 0
        for records in transcript.itervalues():
            exon_record = records["exon"]
            exon_len += exon_record["end"] - exon_record["start"] + 1

            if protein_coding and ("CDS" in records):
                cds_record = records["CDS"]
                cds_len += cds_record["end"] - cds_record["start"] + 1

        if protein_coding:
            well_formed = (cds_len and (cds_len % 3 == 0))

        if well_formed or options.keep_malformed_proteins:
            selection.append(((well_formed, exon_len), transcript))

    if options.keep_all_transcripts or not selection:
        return [item[-1] for item in selection]

    return [max(selection, key=lambda item: item[0])[-1]]


def _do_build_feature_table(options, table, features, protein_coding):
    def add_records(records):
        for record in records:
            features[record["feature"]].append(record)

    retained = read = 0
    for transcripts in table.itervalues():
        read += len(transcripts)
        for exons in select_transcripts(options, transcripts, protein_coding):
            retained += 1
            add_records(get_introns(exons))
            yield (exons, add_records)

    print("\t- Processed %i transcripts, filtered %i (%.1f%%) ..."
          % (read, read - retained, (100.0 * (read - retained)) / read))


def build_coding_seqs_table(options, table):
    """Takes a table generated from a GTF file, and constructs a table for each
       feature, inferring introns and UTRs from the exons and CDSs of the input
       table."""
    print("Building table of features for coding sequences ...")
    features = {"UTR5": [],
                "UTR3": [],
                "CDS": [],
                "intron": []}

    feature_table = _do_build_feature_table(options, table, features, True)
    for (exons, add_records) in feature_table:
        split_exons(exons, add_records)
    return features


def build_noncoding_seqs_table(options, table):
    print("Building table of features for non-coding sequences ...")
    features = {"exon": [],
                "intron": []}

    feature_table = _do_build_feature_table(options, table, features, False)
    for (exons, add_records) in feature_table:
        add_records(record["exon"] for record in exons.itervalues())
    return features


def write_bed(table, target):
    if not table:
        return

    def sort_key(record):
        return (record["contig"], record["start"], record["end"])

    with open(target, "w") as out:
        for record in sorted(table, key=sort_key):
            out.write("%s\t%i\t%i\t%s\t%i\t%s\n" %
                      # As described on http://genome.ucsc.edu/FAQ/FAQformat
                      (record["contig"],      # chrom
                       record["start"],       # chromStart
                       record["end"] + 1,     # chromEnd, past-the-end
                       record["transcript"],  # name
                       0,                     # score
                       record["strand"]))     # strand


###############################################################################
###############################################################################

def parse_arguments(argv):
    prog = "paleomix gtf_to_bed"
    usage = "%s [options] in.gtf out_prefix [in.scaffolds]" % (prog,)

    parser = ArgumentParser(prog=prog, usage=usage)
    parser.add_argument('infile', metavar="INPUT.gtf",
                        help="Input file in GTF format.")
    parser.add_argument('output_prefix', metavar="OUTPUT_PREFIX",
                        help="Prefix of output files.")
    parser.add_argument('scaffolds', metavar="SCAFFOLDS", nargs="?",
                        help="Mapping of scaffolds to contig positions; e.g. "
                             "mapping individual Un* scaffolds onto chrUn.")
    parser.add_argument("--keep-all-transcripts",
                        action="store_true", default=False,
                        help="Include all transcripts in the output BED "
                             "files, not just the longest transcript of each "
                             "gene [default: off]")
    parser.add_argument("--keep-malformed-proteins",
                        action="store_true", default=False,
                        help="Include transcripts of protein-coding in the "
                             "output, even if the the length of the CDS is "
                             "not divisible by 3 [default: off]")
    parser.add_argument('--contig-prefix', default="",
                        help="Add prefix to contig names (e.g. 'chr') "
                             "[default: no prefix].")

    return parser.parse_args(argv)


def main(argv):
    args = parse_arguments(argv)

    scaffolds = {}
    if args.scaffolds:
        print("Reading scaffolds information from %r" % (args.scaffolds,))
        scaffolds = read_scaffolds(args.scaffolds)

    with open_ro(args.infile) as gtf_file:
        print("Reading GTF from %r" % (args.infile,))
        src_table = read_gtf(gtf_file, scaffolds, args.contig_prefix)

    for (source, table) in src_table.iteritems():
        print("Writing tables for '%s' ..." % source)

        if source.startswith("protein"):
            features = build_coding_seqs_table(args, table)
        else:
            features = build_noncoding_seqs_table(args, table)

        for feature in features:
            fpath = "%s.%s.%s.bed" % (args.output_prefix, source, feature)

            print("\tWriting %ss to '%s' ..." % (feature, fpath, ))
            write_bed(features[feature], fpath)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
