.. _bam_pipeline:

BAM Pipeline
============

**Table of Contents:**

.. toctree::

   overview.rst
   requirements.rst
   configuration.rst
   usage.rst
   makefile.rst
   filestructure.rst


The BAM Pipeline is a pipeline designed for the processing of demultiplexed high-throughput sequencing (HTS) data, primarily that generated from Illumina high-throughput sequencing platforms.

The pipeline carries out trimming of adapter sequences, filtering of low quality reads, merging of overlapping mate-pairs to reduce the error rate, mapping of reads using against one or more reference genomes / sequences, filtering of PCR duplicates, analyses of / correction for post-mortem DNA damage, estimation of average coverages and depth-of-coverage histograms, and more. To ensure the correctness of the results, the pipeline invokes frequent validation of intermediate results (*e.g.* using Picard Tools ValidateSamFile.jar), and attempts to detect common errors in input files.

To allow tailoring of the process to the needs of individual projects, many features may be disabled, and the behavior of most programs can be tweaked to suit the specific of a given project, down to and including only carrying out trimming of FASTQ reads, to facilitate use in other contexts.
