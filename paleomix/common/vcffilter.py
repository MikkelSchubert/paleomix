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
from __future__ import with_statement

import sys
import optparse
import collections

import pysam

import paleomix.common.vcfwrap as vcfwrap


_INF = float("inf")
# Rough number of records to keep in memory at once
_CHUNK_SIZE = 10000


def add_varfilter_options(parser):
    group = optparse.OptionGroup(parser, "varFilter: Novel options")
    group.add_option("--homozygous-chromosome", action="append", default=[],
                     help="Filter heterozygous SNPs observed on this "
                          "chromosome (e.g. chrX) %default.")
    group.add_option("-q", "--min-quality", type=int, default=30,
                     help="Minimum Phred score recorded in the QUAL column "
                          "[%default]")
    group.add_option("-f", "--min-allele-frequency", type=float, default=0.2,
                     help="Minimum frequency of the alleles at heterozygous "
                          "sites [%default]. WARNING: A pileup must be "
                          "provided for multi-allelic sites to be filtered!")
    group.add_option("-b", "--pileup", default=None,
                     help="Tabix indexed pileup for multi-allelic sites. This "
                          "is required for such sites to be filtered using "
                          "the --min-allele-frequency filter.")
    group.add_option("-k", "--keep-ambigious-genotypes",
                     default=False, action="store_true",
                     help="Keep SNPs without a most likely genotype "
                          "(based on PL) [%default]")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, "varFilter: Derived options")
    # Options adapted from varFilter
    group.add_option("-Q", "--min-mapping-quality", type=int, default=10,
                     help="Minimum RMS mapping quality for SNPs [%default]")
    group.add_option("-d", "--min-read-depth", type=int, default=8,
                     help="Minimum read depth [%default]")
    group.add_option("-D", "--max-read-depth", type=int, default=10000000,
                     help="Maximum read depth [%default]")
    group.add_option("-a", "--min-num-alt-bases", type=int, default=2,
                     help="Minimum number of alternative bases observed for "
                          "variants [%default]")
    group.add_option("-w", "--min-distance-to-indels", type=int, default=3,
                     help="SNP within INT bp around a gap to be filtered "
                          "[%default]")
    group.add_option("-W", "--min-distance-between-indels",
                     type=int, default=10,
                     help="Window size for filtering adjacent gaps "
                          "[%default]")
    group.add_option("-1", "--min-strand-bias", type=float, default=1e-4,
                     help="Min P-value for strand bias (given PV4) "
                          "[%default]")
    group.add_option("-2", "--min-baseq-bias", type=float, default=1e-100,
                     help="Min P-value for baseQ bias (given PV4) "
                          "[%default]")
    group.add_option("-3", "--min-mapq-bias", type=float, default=0,
                     help="Min P-value for mapQ bias (given PV4) "
                          "[%default]")
    group.add_option("-4", "--min-end-distance-bias", type=float, default=1e-4,
                     help="Min P-value for end distance bias (given PV4) "
                          "[%default]")
    parser.add_option_group(group)


def describe_filters(options):
    return {
        "HET": "Heterozygous SNPs observed on homozygous chromosome (e.g. chrX)",
        "q:%i" % options.min_quality: "Minimum Phred score recorded in the QUAL column",
        "f:%.4f" % options.min_allele_frequency:  "Minimum frequency of the alleles at heterozygous sites",
        "k": "SNPs without a most likely genotype (based on PL)",
        "Q:%i" % options.min_mapping_quality: "Minimum RMS mapping quality",
        "d:%i" % options.min_read_depth: "Minimum read depth",
        "D:%i" % options.max_read_depth: "Maximum read depth",
        "a:%i" % options.min_num_alt_bases: "Minimum number of alternative bases observed for variants",
        "w:%i" % options.min_distance_to_indels: "SNP within INT bp around a gap",
        "W:%i" % options.min_distance_between_indels: "Indel within INT bp of another indel",
        "1:%e" % options.min_strand_bias: "Min P-value for strand bias (given PV4)",
        "2:%e" % options.min_baseq_bias: "Min P-value for baseQ bias (given PV4)",
        "3:%e" % options.min_mapq_bias: "Min P-value for mapQ bias (given PV4)",
        "4:%e" % options.min_end_distance_bias: "Min P-value for end distance bias (given PV4)",
    }


def filter_vcfs(options, vcfs):
    vcfs  = iter(vcfs)
    chunk = collections.deque()
    filename = options.pileup
    min_freq = options.min_allele_frequency

    with AlleleFrequencies(filename, min_freq) as frequencies:
        while _read_chunk(vcfs, chunk):
            chunk = _filter_chunk(options, chunk, frequencies)
            for vcf in _trim_chunk(options, chunk):
                if vcf.filter == ".":
                    vcf.filter = "PASS"

                yield vcf


class AlleleFrequencies:
    VALID, INVALID, NA = range(3)

    def __init__(self, filename, min_freq):
        assert min_freq >= 0
        self._min_freq = min_freq
        self._handle   = None

        if filename and min_freq:
            self._handle = pysam.Tabixfile(filename)
            self.frequency_is_valid = self._frequency_is_valid
        else:
            self.frequency_is_valid = self._frequency_is_always_valid

    def _frequency_is_always_valid(self, contig, position, _ref, _first, _second):
        if self._min_freq:
            sys.stderr.write("WARNING: Multi-allelic SNP found at %s:%i, but --pileup has not been specified.\n" \
                             % (contig, position + 1))
        return self.VALID

    def _frequency_is_valid(self, contig, position, ref, first, second):
        assert self._handle

        if any((len(nt) > 1) for nt in (ref, first, second)):
            assert len(first) != len(second)
            first  = len(first) - len(ref)
            second = len(second) - len(ref)

        counts   = self._fetch(contig, position)
        n_first  = counts.get(first,  0)
        n_second = counts.get(second, 0)

        n_minor = min(n_first, n_second)
        n_major = max(n_first, n_second)
        if not n_major:
            return self.NA
        elif n_minor / float(n_minor + n_major) < self._min_freq:
            return self.INVALID
        return self.VALID

    def close(self):
        if self._handle:
            self._handle.close()
            self._handle = None
        self.frequency_is_valid = self._frequency_is_always_valid

    def __enter__(self):
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback):
        self.close()

    def _fetch(self, contig, position):
        fetched = False
        for line in self._handle.fetch(contig, position, position + 1):
            fields, fetched = line.split("\t"), True
            assert len(fields) == 6
            break

        if not fetched:
            raise RuntimeError("Pileup did not contain position %s:%i, please rebuild." \
                               % (contig, position + 1))
        elif (fields[0] != contig) or (int(fields[1]) != position + 1):
            raise RuntimeError("Got wrong record (%s:%i vs %s:%s), is index corrupt?" \
                               % (contig, position + 1, fields[0], fields[1]))

        counts = {}
        bases  = list(fields[4][::-1].upper())
        ref    = fields[2]
        while bases:
            current = bases.pop()
            if current in "ACGTN":
                counts[current] = counts.get(current, 0) + 1
            elif current in ",.":
                counts[ref] = counts.get(ref, 0) + 1
            elif current in "+-":
                indel_length = [current]
                while bases[-1].isdigit():
                    indel_length.append(bases.pop())
                indel_length = int("".join(indel_length))

                for _ in xrange(abs(indel_length)):
                    bases.pop()

                counts[indel_length] = counts.get(indel_length, 0) + 1
            elif current == "*":
                counts[-1] = counts.get(-1, 0) + 1
            elif current == "^":
                bases.pop()
            elif current != "$":
                raise RuntimeError("Error parsing pileup (unexpected char '%s'): %s" \
                                       % (current, repr(line)))
        return counts


def _read_chunk(vcfs, chunk):
    try:
        while len(chunk) < _CHUNK_SIZE:
            chunk.append(vcfs.next())
    except StopIteration:
        chunk.append(None)

    return len(chunk) > 1


def _trim_chunk(options, chunk):
    min_distance = max(options.min_distance_between_indels,
                       options.min_distance_to_indels)

    if not chunk:
        return
    elif chunk[-1] is None:
        end_chr = "!@#$%^&*()_+"
        end_pos = _INF
        chunk.pop()
    else:
        end_chr = chunk[-1].contig
        end_pos = chunk[-1].pos

    while chunk:
        vcf = chunk[0]
        if  (vcf.contig == end_chr):
            # 'length' will become a too large value for heterozygous SNPs,
            # but it is faster than having to parse every position, and has
            # no effect on the final results.
            length = max(len(vcf.ref), len(vcf.alt))
            if (vcf.pos + length + min_distance) >= end_pos:
                break

        yield chunk.popleft()


def _group_indels_near_position(indels, distance):
    """Returns a dictionary of positions that are either directly covered by, or
    adjacent to indels, given some arbitrary distance. For each position, a list
    of adjacent/overlapping indels are provided."""
    positions = collections.defaultdict(list)
    if not distance:
        return positions

    for vcf in indels:
        # The number of bases covered (excluding the prefix)
        # For ambigious indels (e.g. in low complexity regions), this ensures
        # that the entire region is considered. Note that we do not need to
        # consider the alternative sequence(s)
        length = len(vcf.ref) - 1

        # Inclusive start/end positions for bases that should be blacklisted
        # Note that vcf.pos is the base just before the insertion/deletion
        start = vcf.pos + 1 - distance
        end   = vcf.pos + 1 + distance + length

        for position in xrange(start, end + 1):
            positions[position].append(vcf)

    return positions


def _select_best_indel(indels):
    """Select the highest quality indel, based on the quality,
    prefering low earlier positions above later positions in
    case of ties."""
    def _indel_by_quality_and_position(indel):
        # The negative position is used to select the first
        # of equally quality indels
        return (float(indel.qual), -indel.pos)

    return max(indels, key = _indel_by_quality_and_position)


def _filter_by_indels(options, chunk):
    """Filters a list of SNPs and Indels, such that no SNP is closer to
    an indel than the value set in options.min_distance_to_indels, and
    such that no two indels too close. If two or more indels are within
    this distance, the indel with the highest QUAL score is retained. When
    no unique highest QUAL score exists, an arbitrary indel is retained
    among those indels with the highest QUAL score. SNPs are filtered
    based on prefiltered Indels."""
    indels = [vcf for vcf in chunk if vcfwrap.is_indel(vcf)]

    distance_between = options.min_distance_between_indels
    indel_blacklist  = _group_indels_near_position(indels, distance_between)
    distance_to      = options.min_distance_to_indels
    snp_blacklist    = _group_indels_near_position(indels, distance_to)

    for vcf in chunk:
        if vcfwrap.is_indel(vcf):
            blacklisted = indel_blacklist.get(vcf.pos + 1, [vcf])
            if vcf is not _select_best_indel(blacklisted):
                _mark_as_filtered(vcf, "W:%i" % distance_between)
        elif (vcf.alt != ".") and (vcf.pos in snp_blacklist):
            # TODO: How to handle heterozygous SNPs near
            _mark_as_filtered(vcf, "w:%i" % distance_to)


def _filter_by_properties(options, vcfs, frequencies):
    """Filters a list of SNPs/indels based on the various properties recorded in
    the info column, and others. This mirrors most of the filtering carried out
    by vcfutils.pl varFilter."""
    for vcf in vcfs:
        if float(vcf.qual) < options.min_quality:
            _mark_as_filtered(vcf, "q:%i" % options.min_quality)

        properties = {}
        for field in vcf.info.split(";"):
            if "=" in field:
                key, value = field.split("=")
            else:
                key, value = field, None
            properties[key] = value

        read_depth = float(properties["DP"])
        if options.min_read_depth > read_depth:
            _mark_as_filtered(vcf, "d:%i" % options.min_read_depth)
        elif options.max_read_depth < read_depth:
            _mark_as_filtered(vcf, "D:%i" % options.max_read_depth)

        if "MQ" in properties:
            if float(properties["MQ"]) < options.min_mapping_quality:
                _mark_as_filtered(vcf, "Q:%i" % options.min_mapping_quality)

        if "PV4" in properties:
            pv4 = [float(value) for value in properties["PV4"].split(",")]
            if (pv4[0] < options.min_strand_bias):
                _mark_as_filtered(vcf, "1:%e" % options.min_strand_bias)
            if (pv4[1] < options.min_baseq_bias):
                _mark_as_filtered(vcf, "2:%e" % options.min_baseq_bias)
            if (pv4[2] < options.min_mapq_bias):
                _mark_as_filtered(vcf, "3:%e" % options.min_mapq_bias)
            if (pv4[3] < options.min_end_distance_bias):
                _mark_as_filtered(vcf, "4:%e" % options.min_end_distance_bias)

        if vcf.alt != ".":
            ref_fw, ref_rev, alt_fw, alt_rev = map(int, properties["DP4"].split(","))
            if (alt_fw + alt_rev) < options.min_num_alt_bases:
                _mark_as_filtered(vcf, "a:%i" % options.min_num_alt_bases)

            ml_genotype = vcfwrap.get_ml_genotype(vcf)
            if (ml_genotype == ("N", "N")) and not options.keep_ambigious_genotypes:
                # No most likely genotype
                _mark_as_filtered(vcf, "k")

            if (ml_genotype[0] != ml_genotype[1]):
                if vcf.contig in options.homozygous_chromosome:
                    _mark_as_filtered(vcf, "HET")

                # Filter by frequency of minor allele
                if vcf.ref in ml_genotype:
                    n_minor = min(ref_fw + ref_rev, alt_fw + alt_rev)
                    n_major = max(ref_fw + ref_rev, alt_fw + alt_rev)

                    if (n_minor / float(n_minor + n_major)) < options.min_allele_frequency:
                        _mark_as_filtered(vcf, "f:%.4f" % options.min_allele_frequency)
                else:
                    state = frequencies.frequency_is_valid(vcf.contig, vcf.pos, vcf.ref, *ml_genotype)
                    if state is frequencies.INVALID:
                        _mark_as_filtered(vcf, "f:%.4f" % options.min_allele_frequency)
                    elif state is frequencies.NA:
                        if _mark_as_filtered(vcf, "F:%.4f" % options.min_allele_frequency):
                            sys.stderr.write("WARNING: Could not determine allele-counts for SNP at %s:%s, filtering ...\n" % (vcf.contig, vcf.pos + 1))


def _filter_chunk(options, chunk, frequencies):
    at_end = False
    if chunk[-1] is None:
        at_end = True
        chunk.pop()

    _filter_by_indels(options, chunk)
    _filter_by_properties(options, chunk, frequencies)

    if at_end:
        chunk.append(None)
    return chunk


def _mark_as_filtered(vcf, filter_name):
    if vcf.filter in (".", "PASS"):
        vcf.filter = filter_name
        return True
    elif filter_name not in vcf.filter.split(";"):
        vcf.filter += ";" + filter_name
        return True
