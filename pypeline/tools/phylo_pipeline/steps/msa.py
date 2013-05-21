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

    return FilterSingletonsMetaNode(input_files  = dict((iter(node.output_files).next(), node) for node in afafiles.subnodes),
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
