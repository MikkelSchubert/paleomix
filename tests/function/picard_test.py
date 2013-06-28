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

from pypeline.node import MetaNode
from pypeline.nodes.picard import \
     ValidateBAMNode, \
     BuildSequenceDictNode, \
     MarkDuplicatesNode, \
     MergeSamFilesNode


def _test_validate_bams(config):
    node_params = {"config" : config,
                   "input_bam"     : "tests/data/alignments/library_1.bam",
                   "dependencies"  : config.dependencies}


    standard = ValidateBAMNode(output_log = os.path.join(config.destination, "validate_standard", "log.txt"),
                               **node_params)

    custom   = ValidateBAMNode.customize(output_log = os.path.join(config.destination, "validate_custom", "log.txt"),
                                         **node_params)
    custom.command.set_option("IGNORE_WARNINGS", "True", sep = "=")


    return MetaNode(description  = "ValidateSamFile",
                    dependencies = [standard, custom.build_node()])



def _test_build_sequence_dict(config):
    def _setup(folder):
        root = os.path.join(config.destination, folder)
        fasta = os.path.join(root, "sequence.fasta")
        os.makedirs(root)
        os.symlink(os.path.abspath("tests/data/rCRS.fasta"), fasta)
        return fasta

    node_params = {"config" : config,
                   "dependencies"  : config.dependencies}

    standard_fasta = _setup("seqdict_standard")
    standard = BuildSequenceDictNode(reference = standard_fasta, **node_params)

    custom_fasta = _setup("seqdict_custom")
    custom = BuildSequenceDictNode.customize(reference = custom_fasta, **node_params)
    custom.command.set_option("QUIET", "True", sep = "=")

    return MetaNode(description  = "BuildSequenceDict",
                    dependencies = [standard, custom.build_node()])




def _test_mark_duplicates(config):
    node_params = {"config" : config,
                   "input_bams"    : ("tests/data/alignments/library_1.bam",
                                      "tests/data/alignments/library_2.bam"),
                   "dependencies"  : config.dependencies}


    standard_root = os.path.join(config.destination, "markdup_standard")
    standard = MarkDuplicatesNode(output_bam     = os.path.join(standard_root, "result.bam"),
                                  output_metrics = os.path.join(standard_root, "result.metrics"),
                                  **node_params)

    custom_root = os.path.join(config.destination, "markdup_custom")
    custom   = MarkDuplicatesNode.customize(output_bam     = os.path.join(custom_root, "result.bam"),
                                            output_metrics = os.path.join(custom_root, "result.metrics"),
                                            **node_params)
    custom.command.set_option("VERBOSITY", "DEBUG", sep = "=")


    return MetaNode(description  = "ValidateSamFile",
                    dependencies = [standard, custom.build_node()])



def _test_merge_sam_files(config):
    node_params = {"config" : config,
                   "input_bams"    : ("tests/data/alignments/library_1.bam",
                                      "tests/data/alignments/library_2.bam"),
                   "dependencies"  : config.dependencies}


    output_standard = os.path.join(config.destination, "merge_standard", "result.bam")
    standard = MergeSamFilesNode(output_bam = output_standard,
                                  **node_params)

    output_custom = os.path.join(config.destination, "merge_custom", "result.bam")
    custom   = MergeSamFilesNode.customize(output_bam = output_custom,
                                            **node_params)
    custom.command.set_option("ASSUME_SORTED", "True", sep = "=")


    return MetaNode(description  = "MergeSamFile",
                    dependencies = [standard, custom.build_node()])





def test_picard(config):
    validate = _test_validate_bams(config)
    seqdict  = _test_build_sequence_dict(config)
    rmdup    = _test_mark_duplicates(config)
    merge    = _test_merge_sam_files(config)

    return MetaNode(description = "Picard Tools",
                    dependencies = (validate,
                                    seqdict,
                                    rmdup,
                                    merge))
