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
from __future__ import print_function
from __future__ import with_statement

import optparse
import itertools
import collections

import pysam

import pypeline.common.vcfwrap as vcfwrap
from pypeline.common.samwrap import Tabixfile


_INF = float("inf")
# Rough number of records to keep in memory at once
_CHUNK_SIZE = 10000
# Number of bases to cache when checking mappability
_MAPPABILITY_CACHE  = 500


def add_varfilter_options(parser):
    group = optparse.OptionGroup(parser, "varFilter: Novel options")
    group.add_option("--homozygous-chromosome", action="append", default = [],
                     help = "Filter heterozygous SNPs observed on this chromosome (e.g. chrX) [%default].")
    group.add_option("-q", "--min-quality", type = int, default = 30,
                     help = "Minimum Phred score recorded in the QUAL column [%default]")
    group.add_option("-f", "--max-major-allele-frequency", type = float, default = 0.8,
                     help = "Maximum frequency of the major allele at heterozygous sites [%default]")
    group.add_option("-k", "--keep-ambigious-genotypes", default = False, action = "store_true",
                     help = "Keep SNPs without a most likely genotype (based on PL) [%default]")
    group.add_option("-m", "--filter-by-mappability", default = None,
                     help = "Filter poorly mappable sites using a .bmap (see 'bam_mappability')")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, "varFilter: Derived options")
    # Options adapted from varFilter
    group.add_option("-Q", "--min-mapping-quality", type = int, default = 10,
                     help = "Minimum RMS mapping quality for SNPs [%default]")
    group.add_option("-d", "--min-read-depth", type = int, default = 8,
                     help = "Minimum read depth [%default]")
    group.add_option("-D", "--max-read-depth", type = int, default = 10000000,
                     help = "Maximum read depth [%default]")
    group.add_option("-a", "--min-num-alt-bases", type = int, default = 2,
                     help = "Minimum number of alternative bases observed for variants [%default]")
    group.add_option("-w", "--min-distance-to-indels", type = int, default = 3,
                     help = "SNP within INT bp around a gap to be filtered [%default]")
    group.add_option("-W", "--min-distance-between-indels", type = int, default = 10,
                     help = "Window size for filtering adjacent gaps [%default]")
    group.add_option("-1", "--min-strand-bias", type = float, default = 1e-4,
                     help = "Min P-value for strand bias (given PV4) [%default]")
    group.add_option("-2", "--min-baseq-bias", type = float, default = 1e-100,
                     help = "Min P-value for baseQ bias (given PV4) [%default]")
    group.add_option("-3", "--min-mapq-bias", type = float, default = 0,
                     help = "Min P-value for mapQ bias (given PV4) [%default]")
    group.add_option("-4", "--min-end-distance-bias", type = float, default = 1e-4,
                     help = "Min P-value for end distance bias (given PV4) [%default]")
    parser.add_option_group(group)


def filter_vcfs(options, vcfs):
    vcfs  = iter(vcfs)
    chunk = collections.deque()
    mappability = Mappability(options.filter_by_mappability)

    while _read_chunk(vcfs, chunk):
        chunk = _filter_chunk(options, chunk, mappability)
        for vcf in _trim_chunk(options, chunk):
            if vcf.filter == ".":
                vcf.filter = "PASS"

            yield vcf


class Mappability:
    def __init__(self, filename):
        self._handle = None
        self._contig = ""
        self._start  = -1
        self._end    = -1

        if filename is not None:
            self._handle = pysam.Fastafile(filename)
            self.is_mappable = self._is_mappable


    def is_mappable(self, contig, position):
        return True


    def _is_mappable(self, contig, position):
        assert position # 1-based

        if (contig == self._contig) and (self._start <= position < self._end):
            return (self._cache[position - self._start] == "1")

        self._contig = contig
        self._start  = position
        self._end    = self._start + _MAPPABILITY_CACHE - 1
        self._cache  = self._handle.fetch(contig,
                                          self._start - 1,
                                          self._end)

        return self.is_mappable(contig, position)


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
        end_pos = chunk[-1].pos - min_distance

    while chunk and ((chunk[0].pos < end_pos) or (chunk[0].contig != end_chr)):
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


def _filter_by_indels(options, chunk):
    """Filters a list of SNPs and Indels, such that no SNP is closer to
    an indel than the value set in options.min_distance_to_indels, and
    such that no two indels too close. If two or more indels are within
    this distance, the indel with the highest QUAL score is retained. When
    no unique highest QUAL score exists, an arbitrary indel is retained
    among those indels with the highest QUAL score. SNPs are filtered
    based on prefiltered Indels."""
    indels = set([vcf for vcf in chunk if vcfwrap.is_indel(vcf)])

    distance_between = options.min_distance_between_indels
    indel_blacklist  = _group_indels_near_position(indels, distance_between)
    distance_to      = options.min_distance_to_indels
    snp_blacklist    = _group_indels_near_position(indels, distance_to)

    for vcf in chunk:
        if vcfwrap.is_indel(vcf):
            indels = indel_blacklist.get(vcf.pos + 1, [vcf])
            if vcf != max(indels, key = lambda indel: float(indel.qual)):
                _mark_as_filtered(vcf, "W=%i" % distance_between)
        elif (vcf.alt != ".") and (vcf.pos in snp_blacklist):
            # TODO: How to handle heterozygous SNPs near
            _mark_as_filtered(vcf, "w=%i" % distance_to)


def _filter_by_properties(options, vcfs, mappability):
    """Filters a list of SNPs/indels based on the various properties recorded in
    the info column, and others. This mirrors most of the filtering carried out
    by vcfutils.pl varFilter."""
    for vcf in vcfs:
        if float(vcf.qual) < options.min_quality:
            _mark_as_filtered(vcf, "q=%i" % options.min_quality)

        properties = {}
        for field in vcf.info.split(";"):
            if "=" in field:
                key, value = field.split("=")
            else:
                key, value = field, None
            properties[key] = value

        read_depth = float(properties["DP"])
        if options.min_read_depth > read_depth:
            _mark_as_filtered(vcf, "d=%i" % options.min_read_depth)
        elif options.max_read_depth < read_depth:
            _mark_as_filtered(vcf, "D=%i" % options.max_read_depth)

        if "MQ" in properties:
            if float(properties["MQ"]) < options.min_mapping_quality:
                _mark_as_filtered(vcf, "Q=%i" % options.min_mapping_quality)

        if "PV4" in properties:
            pv4 = [float(value) for value in properties["PV4"].split(",")]
            if (pv4[0] < options.min_strand_bias):
                _mark_as_filtered(vcf, "1=%e" % options.min_strand_bias)
            if (pv4[1] < options.min_baseq_bias):
                _mark_as_filtered(vcf, "2=%e" % options.min_baseq_bias)
            if  (pv4[2] < options.min_mapq_bias):
                _mark_as_filtered(vcf, "3=%e" % options.min_mapq_bias)
            if (pv4[3] < options.min_end_distance_bias):
                _mark_as_filtered(vcf, "4=%e" % options.min_end_distance_bias)

        if vcf.alt != ".":
            # Only filter SNPs?
            if not mappability.is_mappable(vcf.contig, vcf.pos):
                _mark_as_filtered(vcf, "m")

            _, _, alt_fw, alt_rev = properties["DP4"].split(",")
            if (int(alt_fw) + int(alt_rev)) < options.min_num_alt_bases:
                _mark_as_filtered(vcf, "a=%i" % options.min_num_alt_bases)

            ml_genotype = vcfwrap.get_ml_phenotype(vcf)
            if (ml_genotype == ("N", "N")) and not options.keep_ambigious_genotypes:
                # No most likely genotype
                _mark_as_filtered(vcf, "k")

            if (ml_genotype[0] != ml_genotype[1]):
                if vcf.contig in options.homozygous_chromosome:
                    _mark_as_filtered(vcf, "HET")
                if float(properties["AF1"]) > options.max_major_allele_frequency:
                    _mark_as_filtered(vcf, "f=%.4f" % options.max_major_allele_frequency)


def _filter_chunk(options, chunk, mappability):
    at_end = False
    if chunk[-1] is None:
        at_end = True
        chunk.pop()

    _filter_by_indels(options, chunk)
    _filter_by_properties(options, chunk, mappability)

    if at_end:
        chunk.append(None)
    return chunk


def _mark_as_filtered(vcf, filter_name):
    if vcf.filter in (".", "PASS"):
        vcf.filter = filter_name
    elif filter_name not in vcf.filter.split(";"):
        vcf.filter += ";" + filter_name
