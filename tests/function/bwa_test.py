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
import os

from pypeline.node import \
    MetaNode
from pypeline.nodes.bwa import \
    BWAIndexNode, \
    BWANode


def _bwa_index(config):
    return BWAIndexNode(input_file="tests/data/rCRS.fasta",
                        prefix=os.path.join(config.destination, "rCRS"),
                        dependencies=config.dependencies)


def _bwa_se(config, index, algorithm):
    node_params = {"input_file_1": "tests/data/sim_reads/mate_1.fastq.gz",
                   "prefix": os.path.join(config.destination, "rCRS"),
                   "reference": "tests/data/rCRS.fasta",
                   "algorithm": algorithm,
                   "dependencies": (config.dependencies, index)}
    template = os.path.join(config.destination, "%s_se_%s", "output.bam")

    std_bam = template % (algorithm, "standard")
    standard = BWANode(output_file=std_bam, **node_params).build_node()

    cus_bam = template % (algorithm, "custom")
    custom = BWANode(output_file=cus_bam, **node_params)
    custom.commands["convert"].add_option("--rg-id", "myRG")
    custom.commands["convert"].add_option("--rg", "PL:Illumina")
    custom = custom.build_node()

    return MetaNode(description="BWA %s SE" % algorithm,
                    dependencies=[standard, custom])


def _bwa_pe(config, index, algorithm):
    node_params = {"input_file_1": "tests/data/sim_reads/mate_1.fastq.gz",
                   "input_file_2": "tests/data/sim_reads/mate_2.fastq.gz",
                   "prefix": os.path.join(config.destination, "rCRS"),
                   "reference": "tests/data/rCRS.fasta",
                   "algorithm": algorithm,
                   "dependencies": (config.dependencies, index)}
    template = os.path.join(config.destination, "%s_pe_%s", "output.bam")

    std_bam = template % (algorithm, "standard")
    standard = BWANode(output_file=std_bam, **node_params).build_node()

    cus_bam = template % (algorithm, "custom")
    custom = BWANode(output_file=cus_bam, **node_params)
    custom.commands["convert"].add_option("--rg-id", "myRG")
    custom.commands["convert"].add_option("--rg", "PL:Illumina")
    custom = custom.build_node()

    return MetaNode(description="BWA %s PE" % algorithm,
                    dependencies=[standard, custom])


def test_bwa(config):
    index = _bwa_index(config)
    dependencies = []
    for algorithm in ("backtrack", "bwasw", "mem"):
        dependencies.append(_bwa_se(config, index, algorithm))
        dependencies.append(_bwa_pe(config, index, algorithm))

    return MetaNode(description="BWA",
                    dependencies=dependencies)
