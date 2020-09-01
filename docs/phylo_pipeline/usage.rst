.. highlight:: Yaml
.. _phylo_usage:

Pipeline usage
==============

The 'phylo\_pipeline mkfile' command can be used to create a makefile template, as with the 'bam\_pipeline mkfile' command (see section :ref:`bam_usage`). This makefile is used to specify the samples, regions of interest (to be analysed), and options for the various programs:

.. code-block:: bash

    $ paleomix phylo mkfile > makefile.yaml

Note that filenames are not specified explicitly with this pipeline, but are instead inferred from the names of samples, prefixes, etc. as described below.

To execute the pipeline, a command corresponding to the step to be invoked is used (see below):

.. code-block:: bash

    $ paleomix phylo <STEP> [OPTIONS] <MAKEFILE>


Samples
-------

The phylogenetic pipeline expects a number of samples to be specified. Each sample has a name and a sex::

    Samples:
      <GROUP>:
        SAMPLE_NAME:
          Sex: ...

Sex is required, and is used to filter SNPs at homozygous sex chromsomes (e.g. chrX and chrY for male humans). Any names may be used, and can simply be set to e.g. 'NA' in case this feature is not used.

Groups are optional, and may be used either for the sake of the reader, or to specify a group of samples in lists of samples, e.g. when excluding samples form a subsequent step, when filtering singletons, or when rooting phylogenetic trees (see below)

For a given sample with name S, and a prefix with name P, the pipeline will expect files to be located at ./data/samples/*S*.*P*.bam.


Regions of interest
-------------------

Analysis is carried out for a set of "Regions of Interest", which is defined a set of named regions specified using BED files::

    RegionsOfInterest:
      NAME:
        Prefix: NAME_OF_PREFIX
        ProteinCoding: yes/no
        IncludeIndels: yes/no

The options 'ProteinCoding' and 'IncludeIndels' takes values 'yes' and 'no' (without quotation marks), and determines the behavior when calling indels. If 'IncludeIndels' is set to yes, indels are included in the consensus sequence, and if 'ProteinCoding' is set to yes, only indels that are a multiple of 3bp long are included.

The name and the prefix determines the location of the expected BED file and the FASTA file for the prefix: For a region of interest named R, and a prefix named P, the pipeline will expect the BED file to be located at ./data/regions/P.R.bed. The prefix file is expected to be located at ./data/prefixes/P.fasta


Genotyping
----------

Genotyping is done by building a pileup using samtools and calling SNPs / indels using bcftools. The command used for full genotyping is similar to the following command:

.. code-block:: bash

    $ samtools mpileup [OPTIONS] | bcftools view [OPTIONS] -

In addition, SNPs / indels are filtered using the script 'vcf_filter', which is included with the pipeline. This script implements the filteres found in "vcfutils.pl varFilter", with some additions.

Options for either method, including for both "samtools mpileup" and the "bcftools view" command is set using the **Genotyping** section of the makefile, and may be set for all regions of interest (default behavior) or for each set of regions of interest::

    Genotyping:
      Defaults:
        ...

The 'Defaults' key specifies that the options given here apply to all regions of interest; in addition to this key, the name of each set of regions of interest may be used, to set specific values for one set of regions vs. another set. Thus, assuming regions of interest 'ROI\_a' and 'ROI\_b', options may be set as follows::

    Genotyping:
      Defaults:
        ...

      ROI_a:
        ...

      ROI_b:
        ...

For each set of regions of interest named ROI, the final settings are derived by first taking the Defaults, and then overwriting values using the value taken from the ROI section (if one such exists). The following shows how to change values in Defaults for a single ROI::

    Genotyping:
      Defaults:
        --switch: value_a

      ROI_N:
        --switch: value_b

In the above, all ROI except "ROI\_N" will use the switch with 'value\_a', while "ROI\_N" will use 'value\_b'. Executing the 'genotyping' step is described below.

Finally, note the "Padding" option; this option specifies a number of bases to include around each interval in a set of regions of interest. The purpose of this padding is to allow filtering of SNPs based on the distance from indels, in the case where the indels are outside the intervals themselves.


Multiple sequence alignment
---------------------------

Multiple sequence alignment (MSA) is currently carried out using MAFFT, if enabled. Note that it is still nessesary to run the MSA command (see below), even if the multiple sequence alignment itself is disabled (for example in the case where indels are not called in the genotyping step). This is because the MSA step is responsible for generating both the unaligned multi-FASTA files, and the aligned multi-FASTA files. It is nessesary to run the 'genotyping' step prior to running the MSA step (see above).

It is possible to select among the various MAFFT algorithms using the "Algorithm" key, and additionally to specify command-line options for the selected algorithm::

    MultipleSequenceAlignment:
      Defaults:
        Enabled: yes

        MAFFT:
          Algorithm: G-INS-i
          --maxiterate: 1000

Currently supported algorithms are as follows (as described on the `MAFFT website`_):

* mafft - The basic program (mafft)
* auto - Equivalent to command 'mafft --auto'
* fft-ns-1 - Equivalent to the command 'fftns --retree 1'
* fft-ns-2 - Equivalent to the command 'fftns'
* fft-ns-i - Equivalent to the command 'fftnsi'
* nw-ns-i - Equivalent to the command 'nwnsi'
* l-ins-i - Equivalent to the command 'linsi'
* e-ins-i - Equivalent to the command 'einsi'
* g-ins-i - Equivalent to the command 'ginsi'

Command line options are specified as key / value pairs, as shown above for the --maxiterate option, in the same manner that options are specified for the genotyping section. Similarly, options may be specified for all regions of interest ("Defaults"), or using the name of a set of regions of interest, in order to set options for only that set of regions.


Phylogenetic inference
----------------------

Maximum likelyhood Phylogenetic inference is carried out using the ExaML program. A phylogeny consists of a named (subsets of) one or more sets of regions of interest, with individual regions partitioned according to some scheme, and rooted on the midpoint of the tree or one or more taxa::

    PhylogeneticInference:
      PHYLOGENY_NAME:
        ExcludeSamples:
          ...

        RootTreesOn: ...

        PerGeneTrees: yes/no

        RegionsOfInterest:
          REGIONS_NAME:
            Partitions: "111"
            SubsetRegions: SUBSET_NAME

        ExaML:
          Replicates: 1
          Bootstraps: 100
          Model: GAMMA

A phylogeny may exclude any number of samples specified in the Samples region, by listing them under the ExcludeSamples. Furthermore, if groups have been specified for samples (e.g. "<name>"), then these may be used as a short-hand for multiple samples, by using the name of the group including the angle-brackets ("<name>").

Rooting is determined using the RootTreesOn options; if this option is not set, then the resulting trees are rooted on the midpoint of the tree, otherwise it is rooted on the clade containing all the given taxa. If the taxa does not form a monophyletic clade, then rooting is done on the monophyletic clade containing the given taxa.

If PerGeneTrees is set to yes, a tree is generated for every named feature in the regions of interest (e.g. genes), otherwise a super-matrix is created based on all features in all the regions of interest specified for the current phylogeny.

Each phylogeny may include one or more sets of regions of interest, specified under the "RegionsOfInterest", using the same names as those specified under the Project section. Each feature in a set of regions of interest may be partitioned according to position specific scheme. These are specified using a string of numbers (0-9), which is then applied across the selected sequences to determine the model for each position. For example, for the scheme "012" and a given nucleotide sequence, models are applied as follows::

    AAGTAACTTCACCGTTGTGA
    01201201201201201201

Thus, the default partitioning scheme ("111") will use the same model for all positions, and is equivalent to the schemes "1", "11", "1111", etc. Similarly, a per-codon-position scheme may be accomplished using "123" or a similar string. In addition to numbers, the character 'X' may be used to exclude specific positions in an alignment. E.g. to exclude the third position in codons, use a string like "11X". Alternatively, Partitions may be set to 'no' to disable per-feature partitions; instead a single partition is used per set of regions of interest.

The options in the ExaML section specifies the number of bootstrap trees to generate from the original supermatrix, the number of phylogenetic inferences to carry out on the original supermatrix (replicate), and the model used (c.f. the ExaML documentation).

The name (PHYLOGENY_NAME) is used to determine the location of the resulting files, by default ./results/TITLE/phylogenies/NAME/. If per-gene trees are generated, an addition two folders are used, namely the name of the regions of interest, and the name of the gene / feature.

For each phylogeny, the following files are generated:

**alignments.partitions**:

    List of partitions used when running ExaML; the "reduced" file contains the same list of partitions, after empty columns (no called bases) have been excluded.

**alignments.phy**:

    Super-matrix used in conjunction with the list of partitions when calling ExaML; the "reduced" file contains the same matrix, but with empty columns (no bases called) excluded.

**alignments.reduced.binary**:

    The reduced supermatrix / partitions in the binary format used by ExaML.


**bootstraps.newick**:

    List of bootstrap trees in Newick format, rooted as specified in the makefile.


**replicates.newick**:

    List of phylogenies inferred from the full super-matrix, rooted as specified in the makefile.

**replicates.support.newick**:

    List of phylogenies inferred from the full super-matrix, with support values calculated using the bootstrap trees, and rooted as specified in the makefile.


Executing the pipeline
----------------------

The phylogenetic pipeline is excuted similarly to the BAM pipeline, except that a command is provided for each step ('genotyping', 'msa', and 'phylogeny'):

.. code-block:: bash

    $ paleomix phylo <COMMAND> [OPTIONS] <MAKEFILE>

Thus, to execute the genotyping step, the following command is used:

.. code-block:: bash

    $ paleomix phylo genotyping [OPTIONS] <MAKEFILE>

In addition, it is possible to run multiple steps by joining these with the plus-symbol. To run both the 'genotyping' and 'msa' step at the same time, use the following command:

.. code-block:: bash

    $ paleomix phylo genotyping+msa [OPTIONS] <MAKEFILE>


.. _MAFFT website: http://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html