from paleomix.nodes.bwa import (
    BWAAlgorithmNode,
    BWABacktrack,
    BWAIndexNode,
    BWASampe,
    BWASamse,
)


########################################################################################
# Indexing


def test_index_description():
    node = BWAIndexNode(input_file="/path/genome.fasta")

    assert str(node) == "creating BWA index for /path/genome.fasta"


########################################################################################
# BWA mem


def test_bwa_mem_description__se():
    node = BWAAlgorithmNode(
        input_file_1="/path/reads_1.fq.gz",
        output_file="/path/output.bam",
        reference="/path/my_genome.fasta",
    )

    assert str(node) == "aligning '/path/reads_1.fq.gz' onto my_genome using BWA mem"


def test_bwa_mem_description__pe():
    node = BWAAlgorithmNode(
        input_file_1="/path/reads_1.fq.gz",
        input_file_2="/path/reads_2.fq.gz",
        output_file="/path/output.bam",
        reference="/path/my_genome.fasta",
    )

    assert str(node) == "aligning '/path/reads_[12].fq.gz' onto my_genome using BWA mem"
