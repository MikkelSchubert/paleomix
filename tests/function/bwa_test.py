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
from pypeline.nodes.bwa import *


def _bwa_index(config):
    return BWAIndexNode(input_file   = "tests/data/rCRS.fasta",
                        prefix       = os.path.join(config.destination, "rCRS"),
                        dependencies = config.dependencies)


def _bwa_aln_se(config, index):
    node_params = {"input_file"   : "tests/data/sim_reads/mate_1.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "reference"    : "tests/data/rCRS.fasta",
                   "dependencies" : (config.dependencies, index)}
    
    standard = SE_BWANode(output_file = os.path.join(config.destination, "aln_se_standard", "output.bam"),
                          **node_params)
    custom = SE_BWANode.customize(output_file = os.path.join(config.destination, "aln_se_custom", "output.bam"),
                                  **node_params)
    custom.commands["samse"].set_parameter("-r", "@RG\tID:1\tPL:Illumina\tPU:123456\tLB:Library_1\tSM:Sample_1")

    return MetaNode(description  = "BWA aln SE",
                    dependencies = [standard, custom.build_node()])



def _bwa_aln_pe(config, index):
    node_params = {"input_file_1" : "tests/data/sim_reads/mate_1.fastq.gz",
                   "input_file_2" : "tests/data/sim_reads/mate_2.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "reference"    : "tests/data/rCRS.fasta",
                   "dependencies" : (config.dependencies, index)}
    
    standard = PE_BWANode(output_file = os.path.join(config.destination, "aln_pe_standard", "output.bam"),
                          **node_params)
    custom   = PE_BWANode.customize(output_file = os.path.join(config.destination, "aln_pe_custom", "output.bam"),
                                    **node_params)
    custom.commands["sampe"].set_parameter("-r", "@RG\tID:1\tPL:Illumina\tPU:123456\tLB:Library_1\tSM:Sample_1")
    
    return MetaNode(description  = "BWA aln PE",
                    dependencies = [standard, custom.build_node()])




def _bwa_sw_se(config, index):
    node_params = {"input_file_1" : "tests/data/sim_reads/mate_1.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "reference"    : "tests/data/rCRS.fasta",
                   "dependencies" : (config.dependencies, index)}
    
    standard = BWASWNode(output_file = os.path.join(config.destination, "sw_se_standard", "output.bam"),
                         **node_params)
    custom = BWASWNode.customize(output_file = os.path.join(config.destination, "sw_se_custom", "output.bam"),
                                 **node_params)
    custom.commands["aln"].set_parameter("-z", 10)

    return MetaNode(description  = "BWA SW SE",
                dependencies = [standard, custom.build_node()])


def _bwa_sw_pe(config, index):
    node_params = {"input_file_1" : "tests/data/sim_reads/mate_1.fastq.gz",
                   "input_file_2" : "tests/data/sim_reads/mate_2.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "reference"    : "tests/data/rCRS.fasta",
                   "dependencies" : (config.dependencies, index)}
    
    standard = BWASWNode(output_file = os.path.join(config.destination, "sw_pe_standard", "output.bam"),
                         **node_params)
    custom = BWASWNode.customize(output_file = os.path.join(config.destination, "sw_pe_custom", "output.bam"),
                                 **node_params)
    custom.commands["aln"].set_parameter("-z", 7 )

    return MetaNode(description  = "BWA SW PE",
                    dependencies = [standard, custom.build_node()])



def test_bwa(config):
    index  = _bwa_index(config)
    aln_se = _bwa_aln_se(config, index)
    aln_pe = _bwa_aln_pe(config, index)
    sw_se  = _bwa_sw_se(config, index)
    sw_pe  = _bwa_sw_pe(config, index)

    return MetaNode(description  = "BWA",
                    dependencies = (aln_se,
                                    aln_pe,
                                    sw_se,
                                    sw_pe))
