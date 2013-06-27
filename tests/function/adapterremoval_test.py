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
from pypeline.nodes.adapterremoval import *


def _adapterremoval_se(config):
    node_params = {"input_files"   : ("tests/data/raw_reads/se_reads_R1_001.fastq.gz",
                                      "tests/data/raw_reads/se_reads_R1_002.fastq.gz"),
                   "dependencies"  : config.dependencies}

    standard = SE_AdapterRemovalNode(output_prefix = os.path.join(config.destination, "se_standard"),
                                     **node_params)
    custom   = SE_AdapterRemovalNode.customize(output_prefix = os.path.join(config.destination, "se_custom"),
                                               **node_params)
    custom.command.set_option("--minlength", 30)

    return MetaNode(description  = "AdapterRemoval_SE",
                    dependencies = [standard, custom.build_node()])



def _adapterremoval_pe(config):
    node_params = {"input_files_1"  : ("tests/data/raw_reads/pe_reads_R1_001.fastq.gz",
                                       "tests/data/raw_reads/pe_reads_R1_002.fastq.gz"),
                   "input_files_2"  : ("tests/data/raw_reads/pe_reads_R2_001.fastq.gz",
                                       "tests/data/raw_reads/pe_reads_R2_002.fastq.gz"),
                   "dependencies"   : config.dependencies}

    standard = PE_AdapterRemovalNode(output_prefix = os.path.join(config.destination, "pe_standard"),
                                     **node_params)
    custom   = PE_AdapterRemovalNode.customize(output_prefix = os.path.join(config.destination, "pe_custom"),
                                               **node_params)
    custom.command.set_option("--minlength", 30)

    return MetaNode(description  = "AdapterRemoval_PE",
                    dependencies = [standard, custom.build_node()])


def test_adapterremoval(config):
    rm_se = _adapterremoval_se(config)
    rm_pe = _adapterremoval_pe(config)

    return MetaNode(description = "AdapterRemoval",
                    dependencies = (rm_se, rm_pe))
