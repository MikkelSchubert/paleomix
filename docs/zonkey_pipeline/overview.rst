Overview of analytical steps
============================

Briefly, the Zonkey pipeline can run admixture tests on pre-defined species categories (asses, horses, and zebras) to evaluate the ancestry proportions found in the test samples. F1-hybrids are expected to show a balance mixture of two species ancestries, although this balance can deviate from the 50:50 expectation in case limited genetic information is available. This is accomplished using ADMIXTURE [Alexander2009]_.

The zonkey pipeline additionally builds maximum likelihood phylogenetic trees, using RAxML [Stamatakis2006]_ for mitochondrial sequence data and using Treeemix [Pickrell2012]_ for autosomal data. In the latter case, phylogenetic affinities are reconstructed twice: First considering no migration edges and secondly allowing for one migration edge. This allows for fine-grained testing of admixture between the sample and any of the species represented in the reference panel.

In cases where an admixture signal is found, the location of the sample in the mitochondrial tree allows for the identification of the maternal species contributing to the hybrid being examined. For equids, this is essential to distinguish between possible the hybrid forms, such as distinguishing between mules (|female| horse x |male| donkey F1-hybrid) and hinnies (|male| horse x |female| donkey F1-hybrid).

Analyses are presented in HTML reports, one per sample and one summary report when analyzing multiple samples. Figures are generated in both as PNG and PDF format in order to facilitate use in publications (see :ref:`zonkey_filestructure`).


Individual analytical steps
---------------------------

During a typical analyses, the Zonkey pipeline will proceed through the following major analytical steps:


1. Analyzing nuclear alignments:

    1. Input BAMs are indexed using the equivalent of 'samtools index'.

    2. Nucleotides at sites overlapping SNPs in the reference panel are sampled to produce a pseudo-haploid sequence, one in which transitions are included and one in which transitions are excluded, in order to account for the presence of *post-mortem* deamination causing base substitutions. The resulting tables are processed using PLINK to generate the prerequisite files for further analyses.

    3. PCA plots are generated using SmartPCA from the EIGENSOFT suite of tools for both panels of SNPs (including and excluding transitions).

    4. Admixture estimates are carried out using ADMIXTURE, with a partially supervised approach by assigning each sample in the reference panel to one of either two groups (caballine and non-caballine equids) or three groups (asses, horses, and zebras), and processing the SNP panels including and excluding transitions. The input sample is not assigned to a group.

    5. Migration edges are modeled using TreeMix, assuming either 0 or 1 migration edge; analyses is carried out on both the SNP panel including transitions and on the SNP panel excluding transitions.

    6. PNG and PDF figures are generated for each analytical step; in addition, the the per-chromosome coverage of the nuclear genome is plotted.


1. Analyzing mitochondrial alignments:

    1. Input BAMs are indexed using the equivalent of 'samtools index'.

    2. The majority nucleotide at each position in the BAM is determined, and the resulting sequence is added to the mitochondrial reference multiple sequence alignment included in the reference panel.

    3. A maximum likelihood phylogeny is inferred using RAxML, and the resulting tree is drawn, rooted on the midpoint of the phylogeny.


3. Generating reports and summaries

    1. A HTML report is generated for each sample, summarizing the data used and presenting (graphically) the results of each analysis carried out above. All figures are available as PNG and PDF (each figure in the report links to its PDF equivalent).

    2. If multiple samples were processed, a summary of all samples is generated, which presents the major results in an abbreviated form.


.. note::

	While the above shows an ordered list of steps, the pipeline may interleave individual steps during runtime, and may execute multiple steps in parallel in when running in multi-threaded mode (see :ref:`zonkey_usage` for how to run the Zonkey pipeline using multiple threads).

.. |male|   unicode:: U+02642 .. MALE
.. |female|   unicode:: U+02640 .. FEMALE
