.. _phylo_pipeline:

Phylogenetic Pipeline
=====================

**Table of Contents:**

.. toctree::

   overview.rst
   requirements.rst
   configuration.rst
   usage.rst
   makefile.rst
   filestructure.rst


.. warning::

    This section of the documentation is currently undergoing a complete rewrite, and may therefore be incomplete in places.


The Phylogenetic Pipeline is a pipeline designed for processing of (one or more) BAMs, in order to carry out genotyping of a set of regions of interest. Following genotyping, multiple sequence alignment may optionally be carried out (this is required if indels were called), and phylogenetic inference may be done on the regions of interest, using a supermatrix approach through ExaML.

Regions of interest, as defined for the Phylogenetic pipeline, are simply any set of regions in a reference sequence, and may span anything from a few short genomic regions, to the complete exome of complex organisms (tens of thousands of genes), and even entire genomes.

While the Phylogenetic pipeline is designed for ease of use in conjunction with the BAM pipeline, but can be used on arbitrary BAM files, provided that these follow the expected naming scheme (see the :ref:`phylo_usage` section).
