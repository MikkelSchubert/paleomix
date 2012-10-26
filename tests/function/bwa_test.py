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
from pypeline.nodes.bwa import *


def _bwa_index(config):
    return BWAIndexNode(input_file   = "tests/data/rCRS.fasta",
                        prefix       = os.path.join(config.destination, "rCRS"),
                        dependencies = config.dependencies)


def _bwa_aln_se(config, index):
    node_params = {"input_file"   : "tests/data/sim_reads/mate_1.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "read_group"   : "@RG\tID:1\tPL:Illumina\tPU:123456\tLB:Library_1\tSM:Sample_1",
                   "dependencies" : (config.dependencies, index)}
    
    standard = SE_BWANode(output_file = os.path.join(config.destination, "aln_se_standard", "output.bam"),
                          **node_params)
    custom   = Node(description = "Custom SE_BWANode: TODO")

    return Node(description  = "BWA aln SE",
                dependencies = [standard, custom])



def _bwa_aln_pe(config, index):
    node_params = {"input_file_1" : "tests/data/sim_reads/mate_1.fastq.gz",
                   "input_file_2" : "tests/data/sim_reads/mate_2.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "read_group"   : "@RG\tID:1\tPL:Illumina\tPU:123456\tLB:Library_1\tSM:Sample_1",
                   "dependencies" : (config.dependencies, index)}
    
    standard = PE_BWANode(output_file = os.path.join(config.destination, "aln_pe_standard", "output.bam"),
                          **node_params)
    custom   = Node(description = "Custom PE_BWANode: TODO")
    
    return Node(description  = "BWA aln PE",
                dependencies = [standard, custom])




def _bwa_sw_se(config, index):
    node_params = {"input_file_1" : "tests/data/sim_reads/mate_1.fastq.gz",
                   "prefix"       : os.path.join(config.destination, "rCRS"),
                   "dependencies" : (config.dependencies, index)}
    
    standard = BWASWNode(output_file = os.path.join(config.destination, "sw_se_standard", "output.bam"),
                         **node_params)
    custom   = Node(description = "Custom BWASWNode: TODO")

    return Node(description  = "BWA SW SE",
                dependencies = [standard, custom])



def test_bwa(config):
    index  = _bwa_index(config)
    aln_se = _bwa_aln_se(config, index)
    aln_pe = _bwa_aln_pe(config, index)
    sw_se  = _bwa_sw_se(config, index)

    return Node(description  = "BWA",
                dependencies = (aln_se,
                                aln_pe,
                                sw_se))
