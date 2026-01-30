# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
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
