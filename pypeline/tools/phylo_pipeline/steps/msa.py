#!/usr/bin/python

import os
import sys

from pypeline import Pypeline
from pypeline.nodes.sequences import CollectSequencesNode, \
    FilterSingletonsMetaNode
from pypeline.nodes.mafft import MetaMAFFTNode

import common
    


def build_msa_nodes(options, settings, interval, taxa, filtering, dependencies):
    if settings["Default"].lower() != "mafft":
        raise RuntimeError("Only MAFFT support has been implemented!")

    sequencedir = os.path.join(options.destination, "alignments", interval["Name"])
    filenames   = common.collect_fasta_files(options, interval, taxa)
    sequences   = common.collect_sequences(options, interval, taxa)

    fastafiles = CollectSequencesNode(fasta_files  = filenames,
                                      destination  = sequencedir,
                                      sequences    = sequences,
                                      dependencies = dependencies)

    afafiles   = MetaMAFFTNode(rootdir      = sequencedir,
                               sequences    = sequences,
                               preset       = settings["MAFFT"]["Algorithm"],
                               dependencies = fastafiles)

    if not any(filtering.itervalues()):
        return afafiles

    return FilterSingletonsMetaNode(input_files  = dict((node.output_files[0], node) for node in afafiles.subnodes),
                                    destination  = sequencedir + ".filtered",
                                    filter_by    = filtering,
                                    dependencies = afafiles)


def chain(pipeline, options, makefiles):
    destination = options.destination # Move to makefile
    for makefile in makefiles:
        nodes     = []
        settings  = makefile["MSAlignment"]
        intervals = makefile["Project"]["Intervals"]
        filtering = makefile["Project"]["Filter Singletons"]
        taxa      = makefile["Project"]["Taxa"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        for interval in intervals.itervalues():
            nodes.append(build_msa_nodes(options, settings, interval, taxa, filtering, makefile["Nodes"]))
        makefile["Nodes"] = tuple(nodes)
    options.destination = destination
