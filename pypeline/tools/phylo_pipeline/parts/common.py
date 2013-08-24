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
import os
import pysam

from pypeline.common.text import parse_lines


def get_bed_file(options, interval):
    prefix  = get_prefix(interval)
    return os.path.join(options.intervals_root, prefix + ".bed")


def get_fasta_files(options, interval, taxa):
    fastafiles = {}
    for taxon in taxa.itervalues():
        name      = taxon["Name"]
        prefix    = get_prefix(interval)
        fastafile = os.path.join(options.destination, "genotypes", "%s.%s.fasta" % (name, prefix))

        fastafiles[name] = fastafile
    return fastafiles


def get_sequences(options, interval):
    sequences = set()
    with open(get_bed_file(options, interval)) as bedhandle:
        for bed in parse_lines(bedhandle, pysam.asBed()):
            sequences.add(bed.name)
    return frozenset(sequences)


def get_prefix(interval):
    return "%s.%s" % (interval["Genome"], interval["Name"])
