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

from optparse import OptionParser

import pypeline.tools.phylo_pipeline.makefile
import pypeline.ui as ui
from pypeline.common.text import parse_lines


def parse_options(argv, parser = None):
    parser = OptionParser()
    parser.add_option("--run",                default = False, action="store_true")
    parser.add_option("--verbose",            default = False, action="store_true")
    parser.add_option("--expand-nodes",       default = False, action="store_true")
    parser.add_option("--max-threads",        default = 12, type = int)
    parser.add_option("--temp-root",          default = "./temp")
    parser.add_option("--samples-root",       default = "./data/samples")
    parser.add_option("--intervals-root",     default = "./data/intervals")
    parser.add_option("--genomes-root",       default = "./data/genomes")
    parser.add_option("--destination",        default = "./results")

    (options, args) = parser.parse_args(argv)

    makefiles = pypeline.tools.phylo_pipeline.makefile.read_makefiles(args)
    return options, makefiles


def get_genome_for_interval(interval, taxon):
    default = interval.get("Genome")
    genomes = taxon.get("Genomes", {})

    return genomes.get(interval["Name"], default)

    
def collect_bed_files(options, interval, taxa):
    bedfiles = {}
    for taxon in taxa.itervalues():
        name      = taxon["Name"]
        prefix    = get_prefix(interval, taxon)
        bedfile   = os.path.join(options.intervals_root, prefix + ".bed")
        
        bedfiles[name] = bedfile
    return bedfiles


def collect_fasta_files(options, interval, taxa):
    fastafiles = {}
    for taxon in taxa.itervalues():
        name      = taxon["Name"]
        prefix    = get_prefix(interval, taxon)
        fastafile = os.path.join(options.destination, "genotypes", "%s.%s.fasta" % (name, prefix))

        fastafiles[name] = fastafile
    return fastafiles
        

def collect_sequences(options, interval, taxa):
    bedfiles  = collect_bed_files(options, interval, taxa)
    if len(set(bedfiles.itervalues())) > 1:
        raise RuntimeError("Support for combining different intervals files not implemented!")

    # Same set of sequences for all genomes
    sequences = set()
    with open(bedfiles.itervalues().next()) as bedhandle:
        for bed in parse_lines(bedhandle, pysam.asBed()):
            sequences.add(bed.name)
    seqmap = dict(zip(sequences, sequences))

    return dict((name, dict.fromkeys(taxa, name)) for name in seqmap)


        
def get_prefix(interval, taxon = None):
    if not taxon:
        return "{Genome}.{Name}"
    
    genome  = get_genome_for_interval(interval, taxon)
    return "%s.%s" % (genome, interval["Name"])
