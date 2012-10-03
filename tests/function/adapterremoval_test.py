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

import framework 
from pypeline.nodes.adapterremoval import *


def test_adapterremoval_se(config):
    return SE_AdapterRemovalNode(input_files = ("tests/data/raw_reads/se_reads_R1_001.fastq.gz",
                                                "tests/data/raw_reads/se_reads_R1_002.fastq.gz"),
                                 output_prefix = os.path.join(config.destination, "adapterremoval_se"),
                                 dependencies = config.dependencies)

def test_adapterremoval_pe(config):
    return PE_AdapterRemovalNode(input_files_1 = ("tests/data/raw_reads/pe_reads_R1_001.fastq.gz",
                                                  "tests/data/raw_reads/pe_reads_R1_002.fastq.gz"),
                                 input_files_2 = ("tests/data/raw_reads/pe_reads_R2_001.fastq.gz",
                                                  "tests/data/raw_reads/pe_reads_R2_002.fastq.gz"),
                                 output_prefix = os.path.join(config.destination, "adapterremoval_pe"),
                                 dependencies = config.dependencies)

    






if __name__ == '__main__':
    framework.run(test_adapterremoval_se,
                  test_adapterremoval_pe)
