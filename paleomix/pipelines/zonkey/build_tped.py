#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Copyright (c) 2014 Mikkel Schubert <MikkelSch@gmail.com>
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
import collections
import itertools
import os
import random
import sys
import tarfile
from io import TextIOWrapper

import pysam

import paleomix.common.argparse as argparse
import paleomix.common.bamfiles as bamtools
import paleomix.common.fileutils as fileutils
import paleomix.pipelines.zonkey.database as database
from paleomix.common.sampling import reservoir_sampling
from paleomix.common.sequences import NT_CODES

_TRANSITIONS = frozenset((("C", "T"), ("T", "C"), ("G", "A"), ("A", "G")))


def _filter_records(handle, flags=bamtools.EXCLUDED_FLAGS):
    for record in handle:
        if not record.flag & flags:
            yield record


class DownsampledBAM:
    def __init__(self, handle, downsample, included_references):
        self._records = collections.defaultdict(list)

        references = handle.references
        if len(references) != len(included_references):
            raise ValueError(
                "Length of 'included_references' must match the "
                "number of references in BAM file."
            )

        records = _filter_records(handle)
        for record in reservoir_sampling(records, downsample):
            key = references[record.tid]
            self._records[key].append(record)

        self.references = references

        self._records = dict(self._records)
        for value in self._records.values():
            value.sort(key=lambda rec: rec.pos)

    def fetch(self, chrom):
        return self._records.get(chrom, ())


class GenotypeSites:
    def __init__(self, records):
        last_chrom = None
        sites = []

        # chrom, pos, ref, ..., nucleotides
        for chrom, pos, line in records:
            # Convert pos from 1-based to 0-based (same as BAM.pos)
            sites.append((int(pos) - 1, line, []))

            assert last_chrom is None or chrom == last_chrom, (chrom, last_chrom)
            last_chrom = chrom

        sites.sort()
        self._sites = collections.deque(sites)

    def process(self, records, statistics):
        count_used = 0
        count_total = 0
        sites = self._sites
        for record_id, record in enumerate(records):
            count_total += 1

            # TODO: Check sorted
            while sites and sites[0][0] < record.pos:
                yield sites.popleft()

            sequence = record.seq
            sites_iter = iter(sites).__next__
            alignment_iter = iter(record.get_aligned_pairs()).__next__
            alignment_end = record.aend

            try:
                read_used = False
                site_pos, _, nucleotides = sites_iter()
                query_pos, ref_pos = alignment_iter()

                while alignment_end > site_pos:
                    if ref_pos is None or query_pos is None or site_pos > ref_pos:
                        query_pos, ref_pos = alignment_iter()
                    elif site_pos < ref_pos:
                        site_pos, _, nucleotides = sites_iter()
                    else:
                        assert ref_pos == site_pos, (ref_pos, site_pos)
                        nucleotide = sequence[query_pos]
                        if nucleotide != "N":
                            nucleotides.append((record_id, nucleotide))
                            read_used = True

                        query_pos, ref_pos = alignment_iter()
                        site_pos, _, nucleotides = sites_iter()
            except StopIteration:
                if not sites:
                    break
            finally:
                if read_used:
                    count_used += 1

        while sites:
            yield sites.popleft()

        statistics["n_reads"] += count_total
        statistics["n_reads_used"] += count_used


class GenotypeReader:
    def __init__(self, filename):
        self._tar_handle = tarfile.open(filename)
        self._handle = TextIOWrapper(self._tar_handle.extractfile("genotypes.txt"))
        self._header = self._handle.readline().rstrip("\r\n").split("\t")
        self.samples = self._header[-1].split(";")

    def __iter__(self):
        for chrom, records in itertools.groupby(
            self._read_records(), lambda rec: rec[0]
        ):
            sys.stderr.write("Reading contig %r information\n" % (chrom,))
            yield chrom, GenotypeSites(records)

    def _read_records(self):
        for line in self._handle:
            yield line.rstrip().split("\t", 2)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self._handle.close()
        self._tar_handle.close()


def process_record(
    chrom,
    pos,
    line,
    nucleotides,
    statistics,
    records,
    out_incl_ts=sys.stdout,
    out_excl_ts=sys.stdout,
):
    # Filter reads that have already been used
    nucleotides = [
        (rec_id, nuc) for rec_id, nuc in nucleotides if rec_id not in records
    ]

    if not nucleotides:
        return
    elif len(nucleotides) == 1:
        # Avoid unnecessary random() call in 'random.choice'
        record_id, nucleotide = nucleotides[0]
    else:
        record_id, nucleotide = random.choice(nucleotides)

    # Fields are expected to contain at least 2 columns, the first being the
    # reference nucleotide, and the last being the sample genotypes
    _, encoded_genotypes = line.strip().rsplit("\t", 1)

    genotypes = []
    for encoded_nucleotide in encoded_genotypes:
        decoded_nucleotides = NT_CODES.get(encoded_nucleotide, ())
        if len(decoded_nucleotides) == 1:
            genotypes.append(decoded_nucleotides)
            genotypes.append(decoded_nucleotides)
        elif len(decoded_nucleotides) == 2:
            genotypes.extend(decoded_nucleotides)
        else:
            raise ValueError(
                "Invalid nucleotide, not bi-allelic: %r" % (encoded_nucleotide,)
            )

    if nucleotide not in genotypes:
        # Exclude SNPs not observed in the reference panel
        return

    # Convert from 0-based to 1-based
    pos += 1

    # Chromosome, SNP identifier, (dummy) pos in (centi)Morgans, position
    output = [chrom, "chr{}_{}".format(chrom, pos), "0", str(pos)]
    output.extend(genotypes)
    output.append(nucleotide)
    output.append(nucleotide)
    output = " ".join(output)

    statistics["n_sites_incl_ts"] += 1
    out_incl_ts.write("{}\n".format(output))

    if tuple(set(genotypes)) not in _TRANSITIONS:
        statistics["n_sites_excl_ts"] += 1
        out_excl_ts.write("{}\n".format(output))

    records.add(record_id)


def write_tfam(filename, data, samples, bam_sample):
    with open(filename, "w") as handle:
        for key in samples:
            row = data.samples[key]
            sex = {"MALE": 1, "FEMALE": 2, "NA": 0}[row["Sex"].upper()]
            # Family, Individual, Paternal ID, Maternal ID, Sex, Phenotype
            handle.write("{0} {0} 0 0 {1} -9\n".format(key, sex))

        handle.write("{0} {0} 0 0 0 -9\n".format(bam_sample))


def write_summary(args, filename, statistics):
    with open(filename, "w") as handle:
        handle.write("name: %s\n" % (args.name,))
        handle.write("filename: %s\n" % (os.path.abspath(args.bam),))

        for key in ("n_reads", "n_reads_used", "n_sites_incl_ts", "n_sites_excl_ts"):
            handle.write("%s: %s\n" % (key, statistics.get(key, "MISSING")))


def process_bam(args, data, bam_handle, mapping):
    reverse_mapping = dict(zip(mapping.values(), mapping))
    raw_references = bam_handle.references
    references = [reverse_mapping.get(name, name) for name in raw_references]

    if args.downsample:
        sys.stderr.write("Downsampling to at most %i BAM records\n" % (args.downsample))
        bam_handle = DownsampledBAM(bam_handle, args.downsample, references)

    statistics = {
        "n_reads": 0,
        "n_reads_used": 0,
        "n_sites_incl_ts": 0,
        "n_sites_excl_ts": 0,
    }

    fileutils.make_dirs(args.root)

    with open(os.path.join(args.root, "incl_ts.tped"), "w") as output_incl:
        with open(os.path.join(args.root, "excl_ts.tped"), "w") as output_excl:
            with GenotypeReader(args.database) as reader:
                for ref, sites in reader:
                    records = set()
                    raw_ref = raw_references[references.index(ref)]

                    sys.stderr.write("Reading %r from BAM\n" % (raw_ref,))
                    raw_sites = bam_handle.fetch(raw_ref)
                    for pos, line, nucleotides in sites.process(raw_sites, statistics):
                        process_record(
                            ref,
                            pos,
                            line,
                            nucleotides,
                            out_incl_ts=output_incl,
                            out_excl_ts=output_excl,
                            statistics=statistics,
                            records=records,
                        )

                write_summary(
                    args,
                    os.path.join(args.root, "common.summary"),
                    statistics=statistics,
                )
                write_tfam(
                    os.path.join(args.root, "common.tfam"),
                    data,
                    reader.samples,
                    args.name,
                )


def parse_args(argv):
    parser = argparse.ArgumentParser(prog="paleomix zonkey:tped")
    parser.add_argument(
        "root",
        metavar="output_folder",
        help="Output folder in which output files are "
        "to be placed; is created if it does not "
        "already exist.",
    )
    parser.add_argument("database", help="Zonkey database file.")
    parser.add_argument("bam", help="Sorted BAM file.")
    parser.add_argument(
        "--seed",
        type=int,
        help="RNG seed used when downsampling reads; "
        "defaults to using system time as seed.",
    )
    parser.add_argument(
        "--downsample",
        type=int,
        default=0,
        help="Sample N reads from the input BAM file, before building the TPED file. "
        "If not set, or set to zero, all reads are used",
    )
    parser.add_argument(
        "--name", default="Sample", help="Name of sample to be used in output."
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    random.seed(args.seed)

    print("Reading reference information from %r" % (args.database,))

    try:
        data = database.ZonkeyDB(args.database)
    except database.ZonkeyDBError as error:
        sys.stderr.write(
            "Error reading database file %r:\n%s\n" % (args.database, error)
        )
        return 1

    with pysam.AlignmentFile(args.bam) as bam_handle:
        bam_info = data.validate_bam_handle(bam_handle)
        if not bam_info:
            return 1
        elif not bam_info.is_nuclear:
            sys.stderr.write(
                "ERROR: BAM file does not contain identifiable nuclear alignments.\n"
            )
            return 1

        process_bam(args, data, bam_handle, bam_info.nuclear_contigs)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
