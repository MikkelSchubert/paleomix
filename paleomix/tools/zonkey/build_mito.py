#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import argparse
import os
import sys

import pysam

from paleomix.common.formats.fasta import FASTA
from paleomix.common.formats.msa import MSA
from paleomix.common.utilities import fragment
from paleomix.common.formats.phylip import interleaved_phy

import paleomix.ui as ui
import paleomix.tools.zonkey.database as database


def majority_base(site):
    best, count = "N", 0
    for nuc in "ACGT":
        if site[nuc] > count:
            best = nuc
            count = site[nuc]
        elif site[nuc] == count:
            best = "N"

    return best


def majority_sequence(handle, padding, contig_name, contig_length):
    sequence = [dict.fromkeys("ACGTN", 0) for _ in xrange(contig_length)]

    for column in handle.pileup(contig_name):
        position = sequence[column.pos]

        for alignment in column.pileups:
            seq = alignment.alignment.seq
            pos = alignment.query_position

            if pos is not None:
                position[seq[pos]] += 1

    if padding:
        offset = len(sequence) - padding
        for idx in xrange(padding):
            dst = sequence[idx]
            src = sequence[idx + offset]

            for key, value in src.iteritems():
                dst[key] += value

        del sequence[-padding:]

    covered = coverage = 0
    for counts in sequence:
        total = sum(counts.itervalues()) - counts["N"]
        coverage += total

        if total:
            covered += 1

    statistics = {
        "sequence_len": len(sequence),
        "sequence_name": contig_name,
        "nucleotides": coverage,
        "covered_sites": covered,
        "covered_pct": round((100.0 * covered) / len(sequence), 1),
        "mean_coverage": round(coverage / float(len(sequence)), 1),
    }

    return statistics, "".join(map(majority_base, sequence))


def align_majority(reference, majority):
    aligned = []
    reference_iter = iter(reference).next

    for nucleotide in majority:
        reference = reference_iter()
        while reference == "-":
            reference = reference_iter()
            aligned.append("-")

        aligned.append(nucleotide)

    return "".join(aligned)


def truncate_sequences(sequences, name):
    result = {}
    to_len = len(sequences[name].sequence)
    for name, record in sequences.iteritems():
        result[name] = FASTA(name=record.name,
                             meta=record.meta,
                             sequence=record.sequence[:to_len])

    return result


def filter_sequences(sequences):
    selection = {}
    for key, record in sequences.iteritems():
        if record.meta is not None:
            if "EXCLUDE" in map(str.strip, record.meta.upper().split(";")):
                continue

        selection[key] = record

    return selection


def sequences_to_msa(sequences):
    records = []
    for name, record in sorted(sequences.iteritems()):
        records.append(record)

    return MSA(records)


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('bam')
    parser.add_argument('output_prefix')

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    data = database.ZonkeyDB(args.database)
    sequences = data.mitochondria

    try:
        handle = pysam.Samfile(args.bam)
    except (IOError, ValueError), error:
        ui.print_err("Error reading BAM file: %s" % (error,))
        return 1

    with handle:
        bam_info = data.validate_bam_handle(handle)
        if bam_info is None:
            return 1
        elif not bam_info.is_mitochondrial:
            ui.print_err("ERROR: BAM does not contain any known mitochondrial "
                         "sequence found in BAM ..")
            return 1

        reference = sequences[bam_info.mt_contig]
        stats, majority = majority_sequence(handle,
                                            padding=bam_info.mt_padding,
                                            contig_name=bam_info.mt_contig,
                                            contig_length=bam_info.mt_length)

        sequences["Sample"] = FASTA(name="Sample",
                                    meta=None,
                                    sequence=align_majority(reference.sequence,
                                                            majority))

        # Truncate all sequences to match the (now) unpadded sample sequence
        sequences = truncate_sequences(sequences, "Sample")

        sequences = filter_sequences(sequences)

    with open(args.output_prefix + ".summary", "w") as handle:
        stats["filename"] = os.path.abspath(args.bam)

        for key, value in sorted(stats.iteritems()):
            handle.write("{}: {}\n".format(key, value))

    with open(args.output_prefix + ".phy", "w") as handle:
        handle.write(interleaved_phy(sequences_to_msa(sequences)))

    with open(args.output_prefix + ".fasta", "w") as handle:
        for key, record in sorted(sequences.iteritems()):
            handle.write(">{}\n".format(key))
            for line in fragment(60, record.sequence):
                handle.write("{}\n".format(line))

        return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
