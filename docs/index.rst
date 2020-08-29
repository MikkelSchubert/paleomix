
Welcome to PALEOMIX's documentation!
====================================

The PALEOMIX pipelines are a set of pipelines and tools designed to aid the rapid processing of High-Throughput Sequencing (HTS) data: The BAM pipeline processes demultiplexed reads from one or more samples, through sequence processing and alignment, to generate BAM alignment files useful in downstream analyses; the Phylogenetic pipeline carries out genotyping and phylogenetic inference on BAM alignment files, either produced using the BAM pipeline or generated elsewhere; and the Zonkey pipeline carries out a suite of analyses on low coverage equine alignments, in order to detect the presence of F1-hybrids in archaeological assemblages.

The pipelines were originally designed with ancient DNA (aDNA) in mind, and includes several features especially useful for the analyses of ancient samples, but can all be used for the processing of modern samples.

The PALEOMIX pipelines have been published in Nature Protocols; if you make use of PALEOMIX in your work, then please cite

  Schubert M, Ermini L, Sarkissian CD, Jónsson H, Ginolhac A, Schaefer R, Martin MD, Fernández R, Kircher M, McCue M, Willerslev E, and Orlando L. "**Characterization of ancient and modern genomes by SNP detection and phylogenomic and metagenomic analysis using PALEOMIX**". Nat Protoc. 2014 May;9(5):1056-82. doi: `10.1038/nprot.2014.063 <http://dx.doi.org/10.1038/nprot.2014.063>`_. Epub 2014 Apr 10. PubMed PMID: `24722405 <http://www.ncbi.nlm.nih.gov/pubmed/24722405>`_.

The Zonkey pipeline has been published in Journal of Archaeological Science; if you make use of this pipeline in your work, then please cite

  Schubert M, Mashkour M, Gaunitz C, Fages A, Seguin-Orlando A, Sheikhi S, Alfarhan AH, Alquraishi SA, Al-Rasheid KAS, Chuang R, Ermini L, Gamba C, Weinstock J, Vedat O, and Orlando L. "**Zonkey: A simple, accurate and sensitive pipeline to genetically identify equine F1-hybrids in archaeological assemblages**". Journal of Archaeological Science. 2007 Feb; 78:147-157. doi: `10.1016/j.jas.2016.12.005 <http://dx.doi.org/10.1016/j.jas.2016.12.005>`_.

For questions, bug reports, and/or suggestions, please use the PALEOMIX `GitHub tracker <https://github.com/MikkelSchubert/paleomix/issues/>`_.


**Table of Contents:**

.. toctree::
   :maxdepth: 2

   introduction.rst
   installation.rst

   bam_pipeline/index.rst
   phylo_pipeline/index.rst
   zonkey_pipeline/index.rst

   other_tools.rst
   examples.rst

   troubleshooting/index.rst

   yaml.rst
   acknowledgements.rst
   related.rst

   references.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
