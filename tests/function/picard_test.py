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

from pypeline.node import Node
from pypeline.nodes.picard import *


def _test_validate_bams(config):
    node_params = {"config" : config,
                   "input_bam"     : "tests/data/alignments/library_1.bam",
                   "dependencies"  : config.dependencies}


    standard = ValidateBAMNode(output_log = os.path.join(config.destination, "validate_standard", "log.txt"),
                               **node_params)

    custom   = ValidateBAMNode.customize(output_log = os.path.join(config.destination, "validate_custom", "log.txt"),
                                         **node_params)
    custom.command.set_parameter("IGNORE_WARNINGS", "True", sep = "=")


    return Node(description  = "ValidateSamFile", 
                dependencies = [standard, custom.build_node()])

    



def _test_mark_duplicates(config):
    node_params = {"config" : config,
                   "input_bams"    : ("tests/data/alignments/library_1.bam",
                                      "tests/data/alignments/library_2.bam"),
                   "dependencies"  : config.dependencies}


    standard = MarkDuplicatesNode(output_bam     = os.path.join(config.destination, "markdup_standard", "result.bam"),
                                  output_metrics = os.path.join(config.destination, "markdup_standard", "result.metrics"),
                                  **node_params)

    custom   = MarkDuplicatesNode.customize(output_bam     = os.path.join(config.destination, "markdup_custom", "result.bam"),
                                            output_metrics = os.path.join(config.destination, "markdup_custom", "result.metrics"),
                                            **node_params)
    custom.command.set_parameter("VERBOSITY", "DEBUG", sep = "=")


    return Node(description  = "ValidateSamFile", 
                dependencies = [standard, custom.build_node()])



def _test_merge_sam_files(config):
    node_params = {"config" : config,
                   "input_bams"    : ("tests/data/alignments/library_1.bam",
                                      "tests/data/alignments/library_2.bam"),
                   "dependencies"  : config.dependencies}


    standard = MergeSamFilesNode(output_bam     = os.path.join(config.destination, "merge_standard", "result.bam"),
                                  **node_params)

    custom   = MergeSamFilesNode.customize(output_bam     = os.path.join(config.destination, "merge_custom", "result.bam"),
                                            **node_params)
    custom.command.set_parameter("ASSUME_SORTED", "True", sep = "=")


    return Node(description  = "MergeSamFile", 
                dependencies = [standard, custom.build_node()])
    




def test_picard(config):
    validate = _test_validate_bams(config)
    rmdup    = _test_mark_duplicates(config)
    merge    = _test_merge_sam_files(config)
    
    return Node(description = "Picard Tools",
                dependencies = (validate, rmdup, merge))
