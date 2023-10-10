#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

from paleomix.nodes.bowtie2 import Bowtie2IndexNode, Bowtie2Node

########################################################################################
# Indexing


def test_index_description() -> None:
    node = Bowtie2IndexNode(input_file="/path/genome.fasta")

    assert str(node) == "creating Bowtie2 index for /path/genome.fasta"


########################################################################################
# BWA mem


def test_bwa_mem_description__se() -> None:
    node = Bowtie2Node(
        input_file_1="/path/reads_1.fq.gz",
        input_file_2=None,
        output_file="/path/output.bam",
        reference="/path/my_genome.fasta",
    )

    assert str(node) == "aligning '/path/reads_1.fq.gz' onto my_genome using Bowtie2"


def test_bwa_mem_description__pe() -> None:
    node = Bowtie2Node(
        input_file_1="/path/reads_1.fq.gz",
        input_file_2="/path/reads_2.fq.gz",
        output_file="/path/output.bam",
        reference="/path/my_genome.fasta",
    )

    assert str(node) == "aligning '/path/reads_[12].fq.gz' onto my_genome using Bowtie2"
