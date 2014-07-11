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
import copy

from pypeline.common.utilities import \
    get_in, \
    set_in
from pypeline.common.timer import \
    BAMTimer
from pypeline.common.bamfiles import \
    BAMRegionsIter

from pypeline.tools.bam_stats.common import \
    collect_readgroups, \
    collect_references, \
    main_wrapper
from pypeline.tools.bam_stats.coverage import \
    READGROUP_TEMPLATE, \
    write_table


##############################################################################
##############################################################################

def build_region_template(args, handle):
    template = {}
    for key in collect_readgroups(args, handle):
        template[key] = dict(READGROUP_TEMPLATE)
    return template


def get_region_table(counts, region, template):
    subtable = counts.get(region)
    if subtable is not None:
        return subtable

    subtable = copy.deepcopy(template)
    counts[region] = subtable

    return subtable


##############################################################################
##############################################################################

def create_or_get_subtable(table, subtable_key, size):
    subtable = get_in(table, subtable_key)
    if subtable is None:
        subtable = dict(READGROUP_TEMPLATE)
        subtable["Size"] = size
        set_in(table, subtable_key, subtable)
    return subtable


def build_table(args, handle, counts):
    table = {}
    for (key, readgroup) in collect_readgroups(args, handle).iteritems():
        sample = readgroup["SM"]
        library = readgroup["LB"]

        for (reference, size) in collect_references(args, handle).iteritems():
            subtable_key = (args.target_name, sample, library, reference)
            subtable = create_or_get_subtable(table, subtable_key, size)

            statistics = counts.get(reference, {}).get(key, {})
            for (stat, value) in statistics.iteritems():
                subtable[stat] += value

    return table


def print_table(args, handle, counts):
    table = build_table(args, handle, counts)
    write_table(table, args.outfile)


##############################################################################
##############################################################################

def process_record(subtable, record, flags, region):
    qname = record.qname
    if qname.startswith("M_") or qname.startswith("MT_"):
        subtable["Collapsed"] += 1
    elif flags & 0x40:  # first of pair
        subtable["PE_1"] += 1
    elif flags & 0x80:  # second of pair
        subtable["PE_2"] += 1
    else:  # Singleton
        subtable["SE"] += 1

    position = record.pos
    start = region.start
    end = region.end

    for (cigar, num) in record.cigar:
        left = min(max(position, start), end)
        right = min(max(position + num, start), end)
        bases_in_region = right - left
        assert 0 <= bases_in_region <= num

        # 0 = 'M', 1 = 'I', 2 = 'D', 7 = '=', 8 = 'X'
        if cigar in (0, 1, 2, 7, 8):
            if bases_in_region:
                subtable["MID    MM"[cigar]] += bases_in_region

            if cigar != 1:  # Everything but insertions
                position += num
        elif cigar == 3:  # N
            position += num


def process_file(handle, args):
    timer = BAMTimer(handle, step=1000000)

    counts = {}
    region_template = build_region_template(args, handle)
    for region in BAMRegionsIter(handle, args.regions):
        if region.name is None:
            # Trailing unmapped reads
            continue

        name = region.name
        if not args.regions and (handle.nreferences > args.max_contigs):
            name = '<Genome>'

        region_table = get_region_table(counts, name, region_template)
        for (_, records) in region:
            for record in records:
                readgroup = args.get_readgroup_func(record)
                readgroup_table = region_table[readgroup]
                process_record(readgroup_table, record, record.flag, region)
                timer.increment(read=record)
    timer.finalize()

    print_table(args, handle, counts)

    return 0


def main(argv):
    return main_wrapper(process_file, argv, ".coverage")


##############################################################################
##############################################################################

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
