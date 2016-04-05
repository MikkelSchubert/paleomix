.. _introduction:

============
Introduction
============

The PALEOMIX pipeline is a set of pipelines and tools designed to enable the rapid processing of High-Throughput Sequencing (HTS) data from modern and ancient samples. Currently, PALEOMIX consists of 2 major pipelines, and one protocol described in [Schubert2014]_, as well as one as of yet unpublished pipeline:

* **The BAM pipeline** operates on de-multiplexed NGS reads, and carries out the steps necessary to produce high-quality alignments against a reference sequence, ultimately outputting one or more annotated BAM files.

* **The Metagenomic pipeline** is a protocol describing how to carry out metagenomic analyses on reads processed by the BAM pipeline, allowing for the characterisation of the metagenomic population of ancient samples. This protocol makes use of tools included with PALEOMIX.

* **The Phylogenetic pipeline** carries out genotyping, multiple sequence alignment, and phylogenetic inference on a set of regions derived from BAM files (e.g. produced using the BAM Pipeline).

* **The Zonkey Pipeline** is a smaller, experimental pipeline, for the detection of F1 hybrids in equids, based on low coverage nuclear genomes (as few as thousands of aligned reads) and (optionally) mitochondrial DNA.

All pipelines operate through a mix of standard bioinformatics tools, including SAMTools [Li2009b]_, BWA [Li2009a]_, and more, as well as custom scripts written to support the pipelines. The automated pipelines have been designed to run analytical in parallel steps where possible, and to run with minimal user-intervention. To guard against failed steps and to allow easy debugging of failures, all analyses are run in individual temporary folders, all output is logged (though only retained if the command fails), and results are only moved into the destination directory upon successful completion of the given task.

In order to faciliate automatic execution, and to ensure that analyses are documented and can be replicated easily, the BAM and the Phylogenetic Pipelines make use of configuration files (hence-forth "makefiles") in `YAML`_ format ; these are text files which describe a project in terms of input files, settings for programs run as part of the pipeline, and which steps to run. For an overview of the YAML format, refer to the included introduction to :ref:`yaml_intro`, or to the official `YAML`_ website. For a thorough discussion of the makefiles used by either pipeline, please refer to the respective sections of the documentation (*i.e.* :ref:`BAM <bam_makefile>` and :ref:`Phylogentic <phylo_makefile>` pipeline).

.. _YAML: http://www.yaml.org
