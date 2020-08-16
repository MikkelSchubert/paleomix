.. _zonkey_panel:

Reference Panel
===============

The :ref:`zonkey_pipeline` operates using a reference panel of SNPs generated from a selection of extant equid species, including the domestic horse (Equus caballus) and the Przewalski’s wild horse (Equus ferus przewalski); within African asses, the domestic donkey (Equus asinus) and the Somali wild ass (Equus africanus); within Asian asses, the onager (Equus hemionus) and the Tibetan kiang (Equus kiang), and; within zebras: the plains zebra (Equus quagga), the mountains zebra (Equus hartmannae) and the Grevyi zebra (Equus grevyi). These samples were obtained from [Orlando2013]_, [DerSarkissian2015]_, and in particular from [Jonsson2014]_, which published genomes of every remaining extant equid species.

The reference panel has been generated using alignments against the Equus caballus reference nuclear genome (equCab2, via `UCSC`_) and the horse mitochondrial genome (NC\_001640.1, via `NCBI`_). The exact samples used to create the latest version of the reference panel are described below.


Obtaining the reference panel
-------------------------------

The latest version of the Zonkey reference panel (dated 2016-11-01) may be downloaded via the following website:

https://github.com/MikkelSchubert/zonkey/releases/

Once this reference panel has been downloaded, it is strongly recommended that you decompress it using the 'bunzip2' command, since this speeds up several analytical steps (at the cost of about 600 MB of additional disk usage). To decompress the reference panel, simply run 'bunzip2' on the file, as shown here:

.. code-block:: bash

    $ bunzip2 database.tar.bz2

.. warning:
    Do not untar the reference panel. The Zonkey pipeline currently expects data files to be stored in a tar archive, and will not work if files have been extracted into a folder. This may change in the future.

Once this has been done, the Zonkey pipeline may be used as described in the :ref:`zonkey_usage` section.


Samples used in the reference panel
-----------------------------------

The following samples have been used in the construction of the latest version of the reference panel:

======  ===================  ======  ===========  =============================
Group   Species              Sex     Sample Name  Publication
======  ===================  ======  ===========  =============================
Horses  *E. caballus*        Male    FM1798       doi:`10.1016/j.cub.2015.08.032 <https://doi.org/10.1016/j.cub.2015.08.032>`_
.       *E. przewalskii*     Male    SB281        doi:`10.1016/j.cub.2015.08.032 <https://doi.org/10.1016/j.cub.2015.08.032>`_
Asses   *E. a. asinus*       Male    Willy        doi:`10.1038/nature12323       <https://doi.org/10.1038/nature12323>`_
.       *E. kiang*           Female  KIA          doi:`10.1073/pnas.1412627111   <https://doi.org/10.1073/pnas.1412627111>`_
.       *E. h. onager*       Male    ONA          doi:`10.1073/pnas.1412627111   <https://doi.org/10.1073/pnas.1412627111>`_
.       *E. a. somaliensis*  Female  SOM          doi:`10.1073/pnas.1412627111   <https://doi.org/10.1073/pnas.1412627111>`_
Zebras  *E. q. boehmi*       Female  BOE          doi:`10.1073/pnas.1412627111   <https://doi.org/10.1073/pnas.1412627111>`_
.       *E. grevyi*          Female  GRE          doi:`10.1073/pnas.1412627111   <https://doi.org/10.1073/pnas.1412627111>`_
.       *E. z. hartmannae*   Female  HAR          doi:`10.1073/pnas.1412627111   <https://doi.org/10.1073/pnas.1412627111>`_
======  ===================  ======  ===========  =============================


Constructing a reference panel
==============================

The following section describes the format used for the reference panel in Zonkey. It is intended for people who are interested in constructing their own reference panels for a set of species.

.. warning::
    At the time of writing, the number of ancestral groups is hardcoded to 2 and 3 groups; support for any number of ancestral groups is planned. Contact me if this is something you need, and I'll prioritize adding this to the Zonkey pipeline.


It is important to note that a reference panel will is created relative to a single reference genome. For example, for the equine reference panel, all alignments and positions are listed relative to the EquCab2.0 reference genome.

The reference consists of a number of files, which are described below:


settings.yaml
-------------

The settings file is a simple YAML-markup file, which species global options that apply to the reference panel. The current setting file looks as follows:

.. code-block:: yaml

    # Database format; is incremented when the format changes
    Format: 1
    # Revision number; is incremented when the database (but not format) changes
    Revision: 20161101
    # Arguments passed to plink
    Plink: "--horse"
    # Number of chromosomes; required for e.g. PCA analyses
    NChroms: 31
    # N bases of padding used for mitochondrial sequences; the last N bases are
    # expected to be the same as the first N bases, in order to allow alignments
    # at this region of the genome, and are combined to generate final consensus.
    MitoPadding: 31
    # The minimum distance between SNPs, assuming an even distribution of SNPs
    # across the genome. Used when --treemix-k is set to 'auto', which is the
    # default behavior. Value from McCue 2012 (doi:10.1371/journal.pgen.1002451).
    SNPDistance: 150000


The *Format* option defines the panel format, reflects the version of the Zonkey pipeline that supports this panel. It should therefore not be changed unless the format, as described on this page, is changed. The *Revision* reflects the version of a specific reference panel, and should be updated every time data or settings in the reference panel is changed. The equid reference panel simply uses the date at which a given version was created as the revision number.

The *Plink* option lists specific options passed to plink. In the above, this includes just the '--horse' option, which specifies the expected number of chromosomes expected for the horse genome and data aligned against the horse genome.

The *NChroms* option specifies the number of autosomal chromosomes for the reference genome used to construct the reference panel. This is requried for running PCA, but will likely be removed in the future (it is redundant due to contigs.txt).

The *MitoPadding* option is used for the mitochondrial reference sequences, and specifies that some number of the bases at the end of the sequences are identical to the first bases in the sequence. Such duplication (or padding) is used to enable alignments spanning the break introduced when representing a circular genome as a FASTA sequence. If no such padding has been used, then this may simply be set to 0.

The *SNPDistance* option is used to calculate the number of SNPs per block when the --treemix-k option is set to 'auto' (the default behavior). This option assumes that SNPs are evenly distributed across the genome, and calculates block size based on the number of SNPs covered for a given sample.


contigs.txt
-----------

The 'contigs.txt' file contains a table describing the chromsomes included in the zonkey analyses:

.. code-block:: text

    ID  Size       Checksum  Ns
    1   185838109  NA        2276254
    2   120857687  NA        1900145
    3   119479920  NA        1375010
    4   108569075  NA        1172002
    5   99680356   NA        1937819
    X   124114077  NA        2499591

The *ID* column specifies the name of the chromosome. Note that these names are expected to be either numerical (i.e. 1, 2, 21, 31) or sex chromosomes (X or Y). The *Size* column must correspond to the length of the chromosome in the reference genome. The *Ns* column, on the other hand, allows for the number of uncalled bases in the reference to be specified. This value is subtracted from the chromosome size when calculating the relative coverage for sex determination.

The *Checksum* column should contain the MD5 sum calculated for the reference sequence or 'NA' if not available. If specified, this value is intended to be compared with the MD5 sums listed in the headers of BAM files analyzed by the Zonkey pipeline, to ensure that the correct reference sequence is used.

.. note::
    This checksum check is currently not supported, but will be added soon.


.. note::
    The mitochondria is not included in this table; only list autosomes to be analyzed.


samples.txt
-----------

The 'samples.txt' table should contains a list of all samples included in the reference panel, and provides various information about these, most important of which is what ancestral groups a given sample belongs to:

.. code-block:: text

    ID    Group(3)  Group(2)      Species         Sex     SampleID    Publication
    ZBoe  Zebra     NonCaballine  E. q. boehmi    Female  BOE         doi:10.1073/pnas.1412627111
    AOna  Ass       NonCaballine  E. h. onager    Male    ONA         doi:10.1073/pnas.1412627111
    HPrz  Horse     Caballine     E. przewalskii  Male    SB281       doi:10.1016/j.cub.2015.08.032


The *ID* column is used as the name of the sample in the text, tables, and figures generated when running the Zonkey pipeline. It is adviced to keep this name short and preferably descriptive about the group to which the sample belongs.

The *Group(2)* and *Group(3)* columns specify the ancestral groups to which the sample belongs, when connsidering either 2 or 3 ancestral groups. Note that Zonkey currently only supports 2 and 3 ancestral groups (see above).

The *Species*, *Sex*, *SampleID*, and *Publication* columns are meant to contain extra information about the samples, used in the reports generated by the Zonkey pipeline, and are not used directly by the pipeline.


mitochondria.fasta
------------------

The 'mitochondria.fasta' file is expected to contain a multi-sequence alignment involving two different set of sequences. Firstly, it must contain one or more reference sequences against which the input mitochondria alignments have been carried out. In addition, it should contain at least one sequence per species in the reference panel.

Zonkey will compare the reference sequences (either or not subtracting the amount of padding specified in the 'settings.txt' file) against the contigs in the input BAM in order to identify mitochondrial sequences. The Zonkey pipeline then uses the alignment of the reference sequence identified to place the sample into the multi-sequence alignment.

By default, all sequences in the 'mitochondria.fasta' file are included in the mitochondrial phylogeny. However, reference sequences can be excluded by adding a 'EXCLUDE' label after the sequence name:

.. code-block:: text

    >5835107Eq_mito3 EXCLUDE
    gttaatgtagcttaataatat-aaagcaaggcactgaaaatgcctagatgagtattctta

Sequences thus marked are not used for the phylogenetic inference itself.


simulations.txt
---------------

The 'simulations.txt' file contains the results of analyzing simulated data sets in order to generate an emperical distribution of deviations from the expected admixture values.

.. code-block:: text

    NReads  K       Sample1         Sample2         HasTS   Percentile   Value
    1000    2       Caballine       NonCaballine    FALSE   0.000        7.000000e-06
    1000    2       Caballine       NonCaballine    FALSE   0.001        1.973480e-04
    1000    2       Caballine       NonCaballine    FALSE   0.002        2.683880e-04
    1000    2       Caballine       NonCaballine    FALSE   0.003        3.759840e-04
    1000    2       Caballine       NonCaballine    FALSE   0.004        4.595720e-04
    1000    2       Caballine       NonCaballine    FALSE   0.005        5.518900e-04
    1000    2       Caballine       NonCaballine    FALSE   0.006        6.591180e-04

The *NReads* column specifies the number of sequence alignments used in the simulated sample (e.g. 1000, 10000, 100000, and 1000000). Zonkey will use these simulations for different numbers of reads to establish lower and upper bounds on the empirical p-values, where the lower bound is selected as the NReads <= to the number of reads analyzed, and the upper bound is selected as the NReads >= to the number of reads analyzed, when running Zonkey.

The *K* column lists the number of ancestral groups specified when the sample was analyzed; in the equine reference panel, this is either 2 or 3.

The *Sample1* and *Sample2* columns lists the two ancestral groups from which the synthetic hybrid was produced. The order in which these are listed does not matter.

The *HasTS* column specifies if transitions were included (TRUE) or excluded (FALSE).

The *Percentile* column specifies the percent of simulations with a *Value* less than or equal to the current *Value*.

The *Value* column lists the absolute observed deviation from the expected admixture proportion (i.e. 0.5).


There is currently no way to generate this automatically table, but having some support for doing this is planned. Note also that zonkey can be run using a hidden option '--admixture-only', which skips all analyses but those required in order to run ADMIXTURE on the data, and thereby makes running ADMIXTURE exactly as it would be run by Zonkey trivial. For example:

   $ paleomix zonkey run --admixture-only database.tar simulation.bam


genotypes.txt
-------------

The 'genotypes.txt' file contains a table of heterozyous sites relative to the reference sequence used for the reference panel.

.. warning::
    Columns in the 'genotypes.txt' file are expected to be in the exact order shown below.


.. code-block:: text

    Chrom  Pos     Ref AAsi;AKia;AOna;ASom;HCab;HPrz;ZBoe;ZGre;ZHar
    1      1094    A   CAACAAAAA
    1      1102    G   AGGAGGGGG
    1      1114    A   AAAAAAAGA
    1      1126    C   CCCCCCCYC
    1      1128    C   CCCCCCCGC
    1      1674    T   GGGGTTGGG
    1      1675    G   GCCGGGGGG


The *Chrom* column is expected to contain only those contigs / chromosomes listed in the 'contigs.txt' file; the *Pos* column contains the 1-based positions of the variable sites relative to the reference sequence. The *Ref* column contains the nucleotide observed in the reference sequence for the current position; it is currently not used, and may be removed in future versions of Zonkey. The final column contains the nucleotides observed for every sample named in 'samples.txt', joined by semi-colons, and a single letter nucleotide for each of these encoded using UIPAC codes (i.e. A equals AA, W equals AT). The equine reference panel does not include sites not called in every sample, but including such sites is possible by setting the nucleotide to 'N' for the sample with missing data.


Packaging the files
-------------------

The reference panel is distributed as a tar archive. For best performance, the files should be laid out so that the genotypes.txt file is the last file in the archive. This may be accomplished with the following command:

.. code-block:: bash

    $ tar cvf database.tar settings.yaml contigs.txt samples.txt mitochondria.fasta simulations.txt examples genotypes.txt

The tar file may be compressed for distribution (bzip2 or gzip), but should be used uncompressed for best performance.


.. _NCBI: https://www.ncbi.nlm.nih.gov/nuccore/5835107
.. _UCSC: https://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Horse&db=0
