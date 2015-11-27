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

from paleomix.atomiccmd.builder import \
     apply_options
from paleomix.nodes.sequences import \
     CollectSequencesNode, \
     FilterSingletonsNode
from paleomix.nodes.mafft import \
     MAFFTNode

import paleomix.common.fileutils as fileutils


def build_msa_nodes(options, settings, regions, filtering, dependencies):
    if settings["Program"].lower() != "mafft":
        raise RuntimeError("Only MAFFT support has been implemented!")

    sequencedir = os.path.join(options.destination, "alignments", regions["Name"])
    # Run on full set of sequences
    sequences = regions["Sequences"][None]

    node = CollectSequencesNode(fasta_files=regions["Genotypes"],
                                destination=sequencedir,
                                sequences=sequences,
                                dependencies=dependencies)

    if settings["Enabled"]:
        fasta_files = {}
        algorithm = settings["MAFFT"]["Algorithm"]
        for sequence in sequences:
            input_file = os.path.join(sequencedir, sequence + ".fasta")
            output_file = os.path.join(sequencedir, sequence + ".afa")

            mafft = MAFFTNode.customize(input_file=input_file,
                                        output_file=output_file,
                                        algorithm=algorithm,
                                        dependencies=node)
            apply_options(mafft.command, settings["MAFFT"])
            fasta_files[output_file] = mafft.build_node()
    else:
        fasta_files = dict((filename, node) for filename in node.output_files)

    if not any(filtering.itervalues()):
        return fasta_files.values()

    destination = sequencedir + ".filtered"
    filtering = dict(filtering)
    filtered_nodes = []

    for (filename, node) in fasta_files.iteritems():
        output_filename = fileutils.reroot_path(destination, filename)
        filtered_node = FilterSingletonsNode(input_file=filename,
                                             output_file=output_filename,
                                             filter_by=filtering,
                                             dependencies=node)

        filtered_nodes.append(filtered_node)

    return filtered_nodes


def chain(_pipeline, options, makefiles):
    destination = options.destination  # Move to makefile
    for makefile in makefiles:
        nodes = []
        settings = makefile["MultipleSequenceAlignment"]
        filtering = makefile["Project"]["FilterSingletons"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        for regions in makefile["Project"]["Regions"].itervalues():
            regions_settings = settings[regions["Name"]]
            nodes.extend(build_msa_nodes(options, regions_settings, regions,
                                         filtering, makefile["Nodes"]))

        makefile["Nodes"] = tuple(nodes)
    options.destination = destination
