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
import argparse
import collections

import paleomix.common.vcfwrap as vcfwrap


_INF = float("inf")
# Rough number of records to keep in memory at once
_CHUNK_SIZE = 10000


def add_varfilter_options(parser):
    parser.add_argument(
        "--homozygous-chromosome",
        action="append",
        default=[],
        help="Filter heterozygous SNPs observed on this chromosome (e.g. chrX)",
    )
    parser.add_argument(
        "-q",
        "--min-quality",
        type=int,
        default=30,
        help="Minimum Phred score recorded in the QUAL column",
    )
    # No longer supported options:
    parser.add_argument("-f", "--min-allele-frequency", help=argparse.SUPPRESS)
    parser.add_argument("-b", "--pileup", help=argparse.SUPPRESS)
    parser.add_argument(
        "-k",
        "--keep-ambigious-genotypes",
        default=False,
        action="store_true",
        help="Keep SNPs without a most likely genotype (based on PL)",
    )

    # Options adapted from varFilter
    parser.add_argument(
        "-Q",
        "--min-mapping-quality",
        type=int,
        default=10,
        help="Minimum RMS mapping quality for SNPs",
    )
    parser.add_argument(
        "-d", "--min-read-depth", type=int, default=8, help="Minimum read depth",
    )
    parser.add_argument(
        "-D", "--max-read-depth", type=int, default=10000000, help="Maximum read depth",
    )
    parser.add_argument(
        "-a",
        "--min-num-alt-bases",
        type=int,
        default=2,
        help="Minimum number of alternative bases observed for variants",
    )
    parser.add_argument(
        "-w",
        "--min-distance-to-indels",
        type=int,
        default=3,
        help="SNP within INT bp around a gap to be filtered",
    )
    parser.add_argument(
        "-W",
        "--min-distance-between-indels",
        type=int,
        default=10,
        help="Window size for filtering adjacent gaps",
    )
    parser.add_argument(
        "-1",
        "--min-strand-bias",
        type=float,
        default=1e-4,
        help="Min P-value for strand bias (given PV4)",
    )
    parser.add_argument(
        "-2",
        "--min-baseq-bias",
        type=float,
        default=1e-100,
        help="Min P-value for baseQ bias (given PV4)",
    )
    parser.add_argument(
        "-3",
        "--min-mapq-bias",
        type=float,
        default=0,
        help="Min P-value for mapQ bias (given PV4)",
    )
    parser.add_argument(
        "-4",
        "--min-end-distance-bias",
        type=float,
        default=1e-4,
        help="Min P-value for end distance bias (given PV4)",
    )
    parser.add_argument_group(parser)


def describe_filters(options):
    return {
        "HET": "Heterozygous SNPs observed on homozygous chromosome (e.g. chrX)",
        "q:%i" % options.min_quality: "Minimum Phred score recorded in the QUAL column",
        "k": "SNPs without a most likely genotype (based on PL)",
        "Q:%i" % options.min_mapping_quality: "Minimum RMS mapping quality",
        "d:%i" % options.min_read_depth: "Minimum read depth",
        "D:%i" % options.max_read_depth: "Maximum read depth",
        "a:%i"
        % options.min_num_alt_bases: "Minimum number of alternative bases observed for variants",
        "w:%i" % options.min_distance_to_indels: "SNP within INT bp around a gap",
        "W:%i"
        % options.min_distance_between_indels: "Indel within INT bp of another indel",
        "1:%e" % options.min_strand_bias: "Min P-value for strand bias (given PV4)",
        "2:%e" % options.min_baseq_bias: "Min P-value for baseQ bias (given PV4)",
        "3:%e" % options.min_mapq_bias: "Min P-value for mapQ bias (given PV4)",
        "4:%e"
        % options.min_end_distance_bias: "Min P-value for end distance bias (given PV4)",
    }


def filter_vcfs(options, vcfs):
    vcfs = iter(vcfs)
    chunk = collections.deque()

    while _read_chunk(vcfs, chunk):
        chunk = _filter_chunk(options, chunk)
        for vcf in _trim_chunk(options, chunk):
            if vcf.filter == ".":
                vcf.filter = "PASS"

            yield vcf


def _read_chunk(vcfs, chunk):
    try:
        while len(chunk) < _CHUNK_SIZE:
            chunk.append(next(vcfs))
    except StopIteration:
        chunk.append(None)

    return len(chunk) > 1


def _trim_chunk(options, chunk):
    min_distance = max(
        options.min_distance_between_indels, options.min_distance_to_indels
    )

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
        if vcf.contig == end_chr:
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
        end = vcf.pos + 1 + distance + length

        for position in range(start, end + 1):
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

    return max(indels, key=_indel_by_quality_and_position)


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
    indel_blacklist = _group_indels_near_position(indels, distance_between)
    distance_to = options.min_distance_to_indels
    snp_blacklist = _group_indels_near_position(indels, distance_to)

    for vcf in chunk:
        if vcfwrap.is_indel(vcf):
            blacklisted = indel_blacklist.get(vcf.pos + 1, [vcf])
            if vcf is not _select_best_indel(blacklisted):
                _mark_as_filtered(vcf, "W:%i" % distance_between)
        elif (vcf.alt != ".") and (vcf.pos in snp_blacklist):
            # TODO: How to handle heterozygous SNPs near
            _mark_as_filtered(vcf, "w:%i" % distance_to)


def _filter_by_properties(options, vcfs):
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

        mapping_qual = properties.get("MQ")
        if mapping_qual is not None and options.min_mapping_quality:
            if mapping_qual == "." or float(mapping_qual) < options.min_mapping_quality:
                _mark_as_filtered(vcf, "Q:%i" % options.min_mapping_quality)

        if "PV4" in properties:
            pv4 = [float(value) for value in properties["PV4"].split(",")]
            if pv4[0] < options.min_strand_bias:
                _mark_as_filtered(vcf, "1:%e" % options.min_strand_bias)
            if pv4[1] < options.min_baseq_bias:
                _mark_as_filtered(vcf, "2:%e" % options.min_baseq_bias)
            if pv4[2] < options.min_mapq_bias:
                _mark_as_filtered(vcf, "3:%e" % options.min_mapq_bias)
            if pv4[3] < options.min_end_distance_bias:
                _mark_as_filtered(vcf, "4:%e" % options.min_end_distance_bias)

        if vcf.alt != ".":
            ref_fw, ref_rev, alt_fw, alt_rev = map(int, properties["DP4"].split(","))
            if (alt_fw + alt_rev) < options.min_num_alt_bases:
                _mark_as_filtered(vcf, "a:%i" % options.min_num_alt_bases)

            ml_genotype = vcfwrap.get_ml_genotype(vcf)
            if (ml_genotype == ("N", "N")) and not options.keep_ambigious_genotypes:
                # No most likely genotype
                _mark_as_filtered(vcf, "k")

            if ml_genotype[0] != ml_genotype[1]:
                if vcf.contig in options.homozygous_chromosome:
                    _mark_as_filtered(vcf, "HET")


def _filter_chunk(options, chunk):
    at_end = False
    if chunk[-1] is None:
        at_end = True
        chunk.pop()

    _filter_by_indels(options, chunk)
    _filter_by_properties(options, chunk)

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
