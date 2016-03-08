Overview of analytical steps
============================

During a typical analyses, the Phylogenetic pipeline will proceed through the following steps.


1. Genotyping

    1. SNPs are called on the provided regions using SAMTools, and the resulting SNPs are filtered using the 'paleomix vcf_filter' tool.

    2. FASTA sequences are constructed from for the regions of interest, using the filtered SNPs generated above, one FASTA file per set of regions and per sample.

2. Multiple sequence alignment

    1. Per-sample files generated in step 1 are collected, and used to build unaligned multi-FASTA files, one per region of interest.

    2. If enabled, multiple-sequence alignment is carried out on these files using MAFFT, to generate aligned multi-FASTA files.

3. Phylogenetic inference

    Following construction of (aligned) multi-FASTA sequences, phylogenetic inference may be carried out using a partioned maximum likelihood approach via ExaML.