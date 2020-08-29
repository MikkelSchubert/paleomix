.. _introduction:

============
Introduction
============

The PALEOMIX pipeline is a set of pipelines and tools designed to enable the rapid processing of High-Throughput Sequencing (HTS) data from modern and ancient samples. Currently, PALEOMIX consists of two pipelines and one protocol described in [Schubert2014]_, as well as one pipeline described in [Schubert2017]_:

* **The BAM pipeline** operates on de-multiplexed NGS reads, and carries out the steps necessary to produce high-quality alignments against a reference sequence, ultimately producing one or more annotated BAM files [Schubert2014]_.

* **The Phylogenetic pipeline** carries out genotyping, multiple sequence alignment, and phylogenetic inference on a set of regions derived from one or more BAM files, such as those BAM files produced using the BAM Pipeline [Schubert2014]_.

* **The Metagenomic protocol** is a protocol describing how to carry out metagenomic analyses on reads processed by the BAM pipeline, allowing for the characterisation of the metagenomic population of ancient samples. This protocol makes use of tools included with PALEOMIX [Schubert2014]_.

* **The Zonkey Pipeline** is a pipeline for the detection of F1 hybrids in equids, based on low coverage nuclear genomes (as few as thousands of aligned reads) and mitochondrial DNA [Schubert2017]_.

All pipelines operate through a mix of standard bioinformatics tools, such as SAMTools [Li2009b]_, BWA [Li2009a]_, in addition to custom scripts written to support the pipelines. The automated pipelines have been designed to run analytical in parallel steps where possible, and to run with minimal user-intervention. To guard against incomplete runs and to allow easy debugging of failures, all analyses are run in individual temporary folders, all output is logged, and results files are only merged into the destination upon successful completion of the given task.

In order to faciliate automatic execution, and to ensure that analyses are documented and can be replicated easily, the BAM and the Phylogenetic Pipelines make use of configuration files (hence-forth "makefiles") in `YAML`_ format. These are text files which describe a project in terms of input files, settings for programs run as part of the pipeline, and which steps to run. For an overview of the YAML format, refer to the included introduction to :ref:`yaml_intro`, or to the official `YAML`_ website. For a thorough discussion of the makefiles used by either pipeline, please refer to the respective sections of the documentation (*i.e.* :ref:`BAM <bam_makefile>` and :ref:`Phylogentic <phylo_makefile>` pipeline).

.. _YAML: http://www.yaml.org
