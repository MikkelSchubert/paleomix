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
from pypeline.nodes.bowtie2 import \
     Bowtie2IndexNode, \
     Bowtie2Node


def _bowtie2_index(config):
    return Bowtie2IndexNode(input_file   = "tests/data/rCRS.fasta",
                            prefix       = os.path.join(config.destination, "rCRS"),
                            dependencies = config.dependencies)


def _bowtie2_aln_se(config, index):
    node_params = {"input_file_1" : "tests/data/sim_reads/mate_1.fastq.gz",
                   "input_file_2" : None,
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "reference"    : "tests/data/rCRS.fasta",
                   "dependencies" : (config.dependencies, index)}

    standard = Bowtie2Node(output_file = os.path.join(config.destination, "aln_se_standard", "output.bam"),
                           **node_params)
    custom   = Bowtie2Node.customize(output_file = os.path.join(config.destination, "aln_se_custom", "output.bam"),
                                     **node_params)
    custom.commands["aln"].set_option("--very-sensitive")

    return MetaNode(description  = "Bowtie2 aln SE",
                    dependencies = [standard, custom.build_node()])



def _bowtie2_aln_pe(config, index):
    node_params = {"input_file_1" : "tests/data/sim_reads/mate_1.fastq.gz",
                   "input_file_2" : "tests/data/sim_reads/mate_2.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "reference"    : "tests/data/rCRS.fasta",
                   "dependencies" : (config.dependencies, index)}

    standard = Bowtie2Node(output_file = os.path.join(config.destination, "aln_pe_standard", "output.bam"),
                           **node_params)
    custom   = Bowtie2Node.customize(output_file = os.path.join(config.destination, "aln_pe_custom", "output.bam"),
                                     **node_params)
    custom.commands["aln"].set_option("--very-sensitive")

    return MetaNode(description  = "Bowtie2 aln PE",
                    dependencies = [standard, custom.build_node()])






def test_bowtie2(config):
    index  = _bowtie2_index(config)
    aln_se = _bowtie2_aln_se(config, index)
    aln_pe = _bowtie2_aln_pe(config, index)

    return MetaNode(description  = "Bowtie2",
                    dependencies = (aln_se,
                                    aln_pe))
