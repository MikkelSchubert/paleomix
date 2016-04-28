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
import sys
import datetime
import itertools
import collections

from paleomix.common.timer import \
    BAMTimer
from paleomix.common.text import \
    padded_table
from paleomix.common.bamfiles import \
    BAMRegionsIter

from paleomix.tools.bam_stats.common import \
    collect_references, \
    collect_readgroups, \
    main_wrapper


##############################################################################
##############################################################################
##

# Maximum depth to record, and hence the number of columns in output table
_MAX_DEPTH = 200
# Maximum number of count patterns (numbers of bases per library for a given
# site) to cache for bulk processing; see MappingsToTotals for implementation
_MAX_CACHE_SIZE = 10000


# Header prepended to output tables
_HEADER = """# Timestamp: %s
#
# Columns:
#   Contig:   Contig, chromosome, or feature for which a depth histogram was
#             created. Unnamed features are named after the chromosome or
#             contig on which they are located, with a star appended. For
#             example "chr1*". If the maximum number of contigs was exceeded,
#             these are collapsed into one meta-contig named "<Genome>".
#   Size:     The total size of the region. Multiple features with the same
#             name are combined into one row, with the size representing the
#             total for these. Note that overlapping bases are counted 2 (or
#             more) times.
#   MaxDepth: Maximum depth to use when calling SNPs, in order to exclude
#             (at least) the 0.5%% most extreme sites based on read depth,
#             not including sites with depth 0.
#   MD_*:     Fraction of sites with a minimum depth of 1-200.
#"""


##############################################################################
##############################################################################

class MappingToTotals(object):
    def __init__(self, totals, region, smlbid_to_smlb):
        self._region = region
        self._totals = totals
        self._map_by_smlbid, self._totals_src_and_dst \
            = self._build_mappings(totals, region.name, smlbid_to_smlb)
        self._cache = collections.defaultdict(int)

    def process_counts(self, counts, last_pos, cur_pos):
        start = self._region.start
        end = self._region.end

        # Pileups tends to contain identical stretches, so
        # try to avoid repeated lookups by aggregating these
        repeats = 1
        last_count = None
        while counts and (last_pos < cur_pos):
            count = counts.popleft()
            if start <= last_pos < end:
                if count == last_count:
                    repeats += 1
                else:
                    if last_count is not None:
                        self._cache[tuple(last_count)] += repeats
                    last_count = count
                    repeats = 1
            last_pos += 1

        if last_count is not None:
            self._cache[tuple(last_count)] += repeats

        if len(self._cache) > _MAX_CACHE_SIZE:
            self.finalize()

    def finalize(self):
        """Process cached counts."""
        for (count, multiplier) in self._cache.iteritems():
            self._update_totals(count, multiplier)
        self._cache.clear()

    def _update_totals(self, count, multiplier=1):
        for (smlbid, count) in enumerate(count):
            if count:
                for lst in self._map_by_smlbid[smlbid]:
                    lst[0] += count

        for (dst_counts, src_count) in self._totals_src_and_dst:
            if src_count[0]:
                dst_counts[src_count[0]] += multiplier
                src_count[0] = 0

    @classmethod
    def _build_mappings(cls, totals, name, smlbid_to_smlb):
        # Accumulators mapped by sample+library IDs
        totals_by_smlbid = [None] * len(smlbid_to_smlb)
        # Accumulators mapped by the corresponding table keys
        totals_by_table_key = {}

        for (smlbid, (sm_key, lb_key)) in enumerate(smlbid_to_smlb):
            keys = [('*', '*', '*'),
                    (sm_key, '*', '*'),
                    (sm_key, '*', name),
                    (sm_key, lb_key, '*'),
                    (sm_key, lb_key, name)]

            mappings = cls._nonoverlapping_mappings(keys, totals,
                                                    totals_by_table_key)
            totals_by_smlbid[smlbid] = mappings

        totals_src_and_dst = []
        for (key, dst) in totals_by_table_key.iteritems():
            totals_src_and_dst.append((totals[key], dst))

        return totals_by_smlbid, totals_src_and_dst

    @classmethod
    def _nonoverlapping_mappings(cls, keys, totals, totals_by_table_key):
        """Returns a tuple of accumulators for a given set of table keys. As
        multiple table keys may share the same accumulator (e.g. if there is
        only one sample, then sample "*" and that sample will be identical),
        the tuple of accumulators may contain fewer items than keys."""

        mapping = []
        totals_used = set()
        for key in keys:
            # Check that accumulator is not already included
            totals_id = id(totals[key])
            if totals_id not in totals_used:
                totals_used.add(totals_id)
                accumulator = totals_by_table_key.setdefault(key, [0])
                mapping.append(accumulator)
        return tuple(mapping)


##############################################################################
##############################################################################

def calc_max_depth(counts):
    counts = dict(counts)
    counts.pop(0, None)

    running_total = sum(counts.values())
    if not running_total:
        return "NA"

    total = float(running_total)
    for (index, count) in sorted(counts.items()):
        # Stop when less than the 0.5% most extreme values are included
        if running_total / total < 0.005:
            # The max is inclusive, so return the depth just before this one
            return index - 1
        running_total -= count

    return "NA"


def print_table(handle, args, totals):
    lengths = collect_references(args, handle)

    if args.outfile == "-":
        output_handle = sys.stdout
    else:
        output_handle = open(args.outfile, "w")

    with output_handle:
        rows = build_table(args.target_name, totals, lengths)
        output_handle.write(_HEADER % datetime.datetime.now().isoformat())
        output_handle.write("\n")
        for line in padded_table(rows):
            output_handle.write(line)
            output_handle.write("\n")


def calculate_depth_pc(counts, length):
    final_counts = [0] * (_MAX_DEPTH + 1)
    for (depth, count) in counts.iteritems():
        final_counts[min(_MAX_DEPTH, depth)] += count

    running_total = sum(final_counts)
    total = float(length)
    for count in final_counts[1:]:
        yield "%.4f" % (running_total / total,)
        running_total -= count


def build_table(name, totals, lengths):
    header = ["Name", "Sample", "Library", "Contig", "Size", "MaxDepth"]
    for index in xrange(1, _MAX_DEPTH + 1):
        header.append("MD_%03i" % (index,))

    rows = [header]
    last_sm = last_lb = None
    for ((sm_key, lb_key, ct_key), counts) in sorted(totals.items()):
        if (sm_key != last_sm) and (last_sm is not None):
            rows.extend("##")
        elif (lb_key != last_lb) and (last_lb is not None):
            rows.append("#")
        last_sm, last_lb = sm_key, lb_key

        if ct_key == "*":
            length = sum(lengths.itervalues())
        else:
            length = lengths[ct_key]

        row = [name, sm_key, lb_key, ct_key, str(length),
               str(calc_max_depth(counts))]
        row.extend(calculate_depth_pc(counts, length))

        rows.append(row)

    return rows


##############################################################################
##############################################################################

def build_key_struct(args, handle):
    structure = collections.defaultdict(set)
    for readgroup in collect_readgroups(args, handle).itervalues():
        lb_key = readgroup["LB"]
        sm_key = readgroup["SM"]
        structure[sm_key].add(lb_key)

    return structure


def build_new_dicts(totals, dst_sm, dst_lb, references):
    totals[(dst_sm, dst_lb, '*')] = collections.defaultdict(int)
    for contig in references:
        totals[(dst_sm, dst_lb, contig)] = collections.defaultdict(int)


def reuse_dicts(totals, dst_sm, dst_lb, src_sm, src_lb, references):
    totals[(dst_sm, dst_lb, '*')] = totals[(src_sm, src_lb, '*')]
    for contig in references:
        totals[(dst_sm, dst_lb, contig)] = totals[(src_sm, src_lb, contig)]


def build_totals_dict(args, handle):
    references = tuple(collect_references(args, handle))
    structure = build_key_struct(args, handle)

    totals = {}
    for (sm_key, libraries) in structure.iteritems():
        for lb_key in libraries:
            if len(references) == 1:
                key = references[0]
                counts = collections.defaultdict(int)
                totals[(sm_key, lb_key, key)] = counts
                totals[(sm_key, lb_key, '*')] = counts
            else:
                build_new_dicts(totals, sm_key, lb_key, references)

        if len(libraries) == 1:
            key = list(libraries)[0]
            reuse_dicts(totals, sm_key, '*', sm_key, key, references)
        else:
            build_new_dicts(totals, sm_key, '*', references)

    if len(structure) == 1:
        key = list(structure)[0]
        reuse_dicts(totals, '*', '*', key, '*', references)
    else:
        build_new_dicts(totals, '*', '*', references)

    return totals


def count_bases(args, counts, record, rg_to_smlbid, template):
    for _ in xrange(record.alen - len(counts)):
        counts.append(list(template))

    key = rg_to_smlbid.get(args.get_readgroup_func(record))
    if key is None:
        # Unknown readgroups are treated as missing readgroups
        key = rg_to_smlbid[None]

    index = 0
    for (cigar, count) in record.cigar:
        if cigar in (0, 7, 8):
            for counter in itertools.islice(counts, index, index + count):
                counter[key] += 1
            index += count
        elif cigar in (2, 3, 6):
            index += count


def build_rg_to_smlbid_keys(args, handle):
    """Returns a dictionary which maps a readgroup ID to an index value,
    as well as a list containing a tuple (samples, library) corresponding
    to each index. Typically, this list will be shorter than the map of read-
    groups, as multiple read-groups will map to the same sample / library.
    """

    rg_to_lbsmid = {}
    lbsm_to_lbsmid = {}
    lbsmid_to_smlb = []
    for (key_rg, readgroup) in collect_readgroups(args, handle).iteritems():
        key_sm = readgroup["SM"]
        key_lb = readgroup["LB"]

        key_lbsm = (key_sm, key_lb)
        if key_lbsm not in lbsm_to_lbsmid:
            lbsm_to_lbsmid[key_lbsm] = len(lbsm_to_lbsmid)
            lbsmid_to_smlb.append(key_lbsm)

        rg_to_lbsmid[key_rg] = lbsm_to_lbsmid[key_lbsm]
    return rg_to_lbsmid, lbsmid_to_smlb


def process_file(handle, args):
    timer = BAMTimer(handle, step=1000000)

    last_tid = 0
    totals = build_totals_dict(args, handle)
    rg_to_smlbid, smlbid_to_smlb = build_rg_to_smlbid_keys(args, handle)
    template = [0] * len(smlbid_to_smlb)

    for region in BAMRegionsIter(handle, args.regions):
        if region.name is None:
            # Trailing unmapped reads
            break
        elif not args.regions and (handle.nreferences > args.max_contigs):
            region.name = '<Genome>'

        last_pos = 0
        counts = collections.deque()
        mapping = MappingToTotals(totals, region, smlbid_to_smlb)
        for (position, records) in region:
            mapping.process_counts(counts, last_pos, position)

            for record in records:
                timer.increment(read=record)
                count_bases(args, counts, record, rg_to_smlbid, template)

            if (region.tid, position) < (last_tid, last_pos):
                sys.stderr.write("ERROR: Input BAM file is unsorted\n")
                return 1

            last_pos = position
            last_tid = region.tid

        # Process columns in region after last read
        mapping.process_counts(counts, last_pos, float("inf"))
        mapping.finalize()
    timer.finalize()

    if not args.ignore_readgroups:
        # Exclude counts for reads with no read-groups, if none such were seen
        for (key, _, _), value in totals.iteritems():
            if key == '<NA>' and value:
                break
        else:
            for key in totals.keys():
                if key[0] == '<NA>':
                    totals.pop(key)

    print_table(handle, args, totals)

    return 0


def main(argv):
    return main_wrapper(process_file, argv, ".depths")

##############################################################################
##############################################################################

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
