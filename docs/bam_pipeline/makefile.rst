.. highlight:: YAML
.. _bam_makefile:

Makefile description
====================

.. contents::

The following sections reviews the options available in the BAM pipeline makefiles. As described in the :ref:`bam_usage` section, a default makefile may be generated using the 'paleomix bam\_pipeline makefile' command. For clarity, the location of options in subsections are specified by concatenating the names using '\:\:' as a separator. Thus, in the following (simplified example), the 'UseSeed' option (line 13) would be referred to as 'Options \:\: Aligners \:\: BWA \:\: UseSeed':

.. code-block:: yaml
    :emphasize-lines: 13
    :linenos:

    Options:

      # Settings for aligners supported by the pipeline
      Aligners:
        # Choice of aligner software to use, either "BWA" or "Bowtie2"
        Program: BWA

        # Settings for mappings performed using BWA
        BWA:
          # May be disabled ("no") for aDNA alignments with the 'aln' algorithm.
          # Post-mortem damage localizes to the seed region, which BWA expects to
          # have few errors (sets "-l"). See http://pmid.us/22574660
          UseSeed: yes


Specifying command-line options
-------------------------------

For several programs it is possible to directly specify command-line options; this is accomplished in one of 3 ways; firstly, for command-line options that take a single value, this is accomplished simply by specifying the option and value as any other option. For example, if we wish to supply the option --mm 5 to AdapterRemoval, then we would list it as "--mm: 5" (all other options omitted for brevity)::

    AdapterRemoval:
      --mm: 5

For options that do not take any values, such as the AdapterRemoval `--trimns` (enabling the trimming of Ns in the reads), these are specified either as "--trimmns: ", with the value left blank, or as "--trimns: yes". The following are therefore equivalent::

    AdapterRemoval:
      --trimns:      # Method 1
      --trimns: yes  # Method 2

In some cases the BAM pipeline will enable features by default, but still allow these to be overridden. In those cases, the feature can be disabled by setting the value to 'no' (without quotes), as shown here::

    AdapterRemoval:
      --trimns: no

If you need to provide the text "yes" or "no" as the value for an option, it is necessary to put these in quotes::

    --my-option: "yes"
    --my-option: "no"

In some cases it is possible or even necessary to specify an option multiple times. Due to the way YAML works, this is not possible to do so directly. Instead, the pipeline allows multiple instances of the same option by providing these as a list::

    --my-option:
        - "yes"
        - "no"
        - "maybe"

The above will be translated into calling the program in question with the options "--my-option yes --my-option no --my-option maybe".


Options
-------

By default, the 'Options' section of the makefile contains the following:

.. literalinclude:: makefile.yaml
    :language: yaml
    :linenos:
    :lines: 2-107


General Options
^^^^^^^^^^^^^^^

**Options \:\: Platform**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 7
        :lines: 7-8

    The sequencing platform used to generate the sequencing data; this information is recorded in the resulting BAM file, and may be used by downstream tools. The `SAM/BAM specification`_ the valid platforms, which currently include 'CAPILLARY', 'HELICOS', 'ILLUMINA', 'IONTORRENT', 'LS454', 'ONT', 'PACBIO', and 'SOLID'.

**Options \:\: QualityOffset**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 9
        :lines: 9-13

    The QualityOffset option refers to the starting ASCII value used to encode `Phred quality-scores`_ in user-provided FASTQ files, with the possible values of 33, 64, and 'Solexa'. For most modern data, this will be 33, corresponding to ASCII characters in the range '!' to 'J'. Older data is often encoded using the offset 64, corresponding to ASCII characters in the range '@' to 'h', and more rarely using Solexa quality-scores, which represent a different scheme than Phred scores, and which occupy the range of ASCII values from ';' to 'h'. For a visual representation of this, refer to the Wikipedia article linked above.


Adapter Trimming
^^^^^^^^^^^^^^^^

The "AdapterRemoval" subsection allows for options that are applied when AdapterRemoval is applied to the FASTQ reads supplied by the user. For a more detailed description of command-line options, please refer to the `AdapterRemoval documentation`_. A few important particularly options are described here:

**Options \:\: AdapterRemoval \:\: --adapter1** and **Options \:\: AdapterRemoval \:\: --adapter2**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 17
        :lines: 17-19


    These two options are used to specify the adapter sequences used to identify and trim reads that contain adapter contamination. Thus, the sequence provided for --adapter1 is expected to be found in the mate 1 reads, and the sequence specified for --adapter2 is expected to be found in the mate 2 reads. In both cases, these should be specified as in the orientation that appear in these files (i.e. it should be possible to grep the files for these, assuming that the reads were long enough, and treating Ns as wildcards). It is very important that these be specified correctly. Please refer to the `AdapterRemoval documentation`_ for more information.


    .. note::
        As of version AdapterRemoval 2.1, it is possible to use multiple threads to speed up trimming of adapter sequences. This is accomplished not by setting the --threads command-line option in the makefile, but by supplying the --adapterremoval-max-threads option to the BAM pipeline itself:

        .. code-block:: bash

            $ paleomix bam run makefile.yaml --adapterremoval-max-threads 2


**Options \:\: AdapterRemoval \:\: --mm**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 22
        :lines: 22

    Sets the fraction of mismatches allowed when aligning reads / adapter sequences. If the specified value (MM) is greater than 1, this is calculated as 1 / MM, otherwise the value is used directly. To set, replace the default value as desired::

        --mm: 3    # Maximum mismatch rate of 1 / 3
        --mm: 5    # Maximum mismatch rate of 1 / 5
        --mm: 0.2  # Maximum mismatch rate of 1 / 5


**Options \:\: AdapterRemoval \:\: --minlength**

    The minimum length required after read merging, adapter trimming, and base-quality quality trimming; resulting reads shorter than this length are discarded, and thereby excluded from further analyses by the pipeline. A value of at least 25 bp is recommended to cut down on the rate of spurious alignments; if possible, a value of 30 bp may be used to greatly reduce the fraction of spurious alignments, with smaller gains for greater minimums [Schubert2012]_.

    .. warning::
        The default value used by PALEOMIX for `--minlength` (25 bp) differs from the default value for AdapterRemoval (15 bp). Thus, if a minimum length of 15 bp is desired, it is nessesarily to explicitly state so in the makefile, simply commenting out this command-line argument is not sufficient.


**Options \:\: AdapterRemoval \:\: --collapse**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 25
        :lines: 25

    If enabled, AdapterRemoval will attempt to combine overlapping paired-end reads into a single (potentially longer) sequence. This has at least two advantages, namely that longer reads allow for less ambiguous alignments against the target reference genome, and that the fidelity of the overlapping region (potentially the entire read) is improved by selecting the highest quality base when discrepancies are observed. The names of reads thus merged are prefixed with either 'M\_' or 'MT\_', with the latter marking reads that have been trimmed from the 5' or 3' termini following collapse, and which therefore do not represent the full insert. To disable this behavior, set the option to 'no' (without quotes)::

        --collapse: yes  # Option enabled
        --collapse: no   # Option disabled

    .. note::
        This option may be combined with the 'ExcludeReads' option (see below), to either eliminate or select for short inserts, depending on the expectations from the experiment. I.e. for ancient samples, where the most inserts should be short enough to allow collapsing (< 2x read read - 11, by default), excluding paired (uncollapsed) and singleton reads may help reduce the fraction of exogenous DNA mapped.


**Options \:\: AdapterRemoval \:\: --trimns**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 26
        :lines: 26

    If set to 'yes' (without quotes), AdapterRemoval will trim uncalled bases ('N') from the 5' and 3' end of the reads. Trimming will stop at the first called base ('A', 'C', 'G', or 'T'). If both --trimns and --trimqualities are enabled, then consecutive stretches of Ns and / or low-quality bases are trimmed from the 5' and 3' end of the reads. To disable, set the option to 'no' (without quotes)::

        --trimns: yes  # Option enabled
        --trimns: no   # Option disabled


**Options \:\: AdapterRemoval \:\: --trimqualities**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 27
        :lines: 27

    If set to 'yes' (without quotes), AdapterRemoval will trim low-quality bases from the 5' and 3' end of the reads. Trimming will stop at the first base which is greater than the (Phred encoded) minimum quality score specified using the command-line option --minquality. This value defaults to 2. If both --trimns and --trimqualities are enabled, then consecutive stretches of Ns and / or low-quality bases are trimmed from the 5' and 3' end of the reads. To disable, set the option to 'no' (without quotes)::

        --trimqualities: yes  # Option enabled
        --trimqualities: no   # Option disabled


Short read aligners
^^^^^^^^^^^^^^^^^^^

This section allow selection between supported short read aligners (currently BWA [Li2009a]_ and Bowtie2 [Langmead2012]_), as well as setting options for these, individually:

.. literalinclude:: makefile.yaml
    :language: yaml
    :linenos:
    :lineno-start: 29
    :lines: 29-32


To select a mapping program, set the 'Program' option appropriately::

    Program: BWA      # Using BWA to map reads
    Program: Bowtie2  # Using Bowtie2 to map reads


Short read aligners - BWA
"""""""""""""""""""""""""

    The following options are applied only when running the BWA short read aligner; see the section "Options: Short read aligners" above for how to select this aligner.

    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 35
        :lines: 35-49

    **Options \:\: Aligners \:\: BWA \:\: Algorithm**
        .. literalinclude:: makefile.yaml
            :language: yaml
            :linenos:
            :lineno-start: 26
            :lines: 36-38

        The mapping algorithm to use; options are 'backtrack' (corresponding to 'bwa aln'), 'bwasw', and 'mem'. Additional command-line options may be specified for these. Algorithms are selected as follows::

            Algorithm: backtrack  # 'Backtrack' algorithm, using the command 'bwa aln'
            Algorithm: bwasw      # 'SW' algorithm for long queries, using the command 'bwa bwasw'
            Algorithm: mem        # 'mem' algorithm, using the command 'bwa mem'

        .. warning::

            Alignment algorithms 'bwasw' and 'mem' currently cannot be used with input data that is encoded using QualityOffset 64 or 'Solexa'. This is a limitation of PALEOMIX, and will be resolved in future versions. In the mean time, this can be circumvented by converting FASTQ reads to the standard quality-offset 33, using for example `seqtk`_.


    **Options \:\: Aligners \:\: BWA \:\: MinQuality**
        .. literalinclude:: makefile.yaml
            :language: yaml
            :linenos:
            :lineno-start: 39
            :lines: 39-40

        Specifies the minimum mapping quality of alignments produced by BWA. Any aligned read with a quality score below this value are removed during the mapping process. Note that while unmapped read have a quality of zero, these are not excluded by a non-zero 'MinQuality' value. To filter unmapped reads, use the option 'FilterUnmappedReads' (see below). To set this option, replace the default value with a desired minimum::

            MinQuality: 0   # Keep all hits
            MinQuality: 25  # Keep only hits where mapping-quality >= 25

    **Options \:\: Aligners \:\: BWA \:\: FilterUnmappedReads**
        .. literalinclude:: makefile.yaml
            :language: yaml
            :linenos:
            :lineno-start: 41
            :lines: 41-42

        Specifies wether or not unmapped reads (reads not aligned to a target sequence) are to be retained in the resulting BAM files. If set to 'yes' (without quotes), all unmapped reads are discarded during the mapping process, while setting the option to 'no' (without quotes) retains these reads in the BAM. By convention, paired reads in which one mate is unmapped are assigned the same chromosome and position, while no chromosome / position are assigned to unmapped single-end reads. To change this setting, replace the value with either 'yes' or 'no' (without quotes)::

            FilterUnmappedReads: yes  # Remove unmapped reads during alignment
            FilterUnmappedReads: no   # Keep unmapped reads

    **Options \:\: Aligners \:\: BWA \:\: \***

        Additional command-line options may be specified for the selected alignment algorithm, as described in the "Specifying command-line options" section above. See also the examples listed for Bowtie2 below. Note that for the 'backtrack' algorithm, it is only possible to specify options for the 'bwa aln' call.



Short read aligners - Bowtie2
"""""""""""""""""""""""""""""
    The following options are applied only when running the Bowtie2 short read aligner; see the section "Options: Short read aligners" above for how to select this aligner.

    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 51
        :lines: 51-65

    **Options \:\: Aligners \:\: Bowtie2 \:\: MinQuality**
        .. literalinclude:: makefile.yaml
            :language: yaml
            :linenos:
            :lineno-start: 53
            :lines: 53-54

        See 'Options \:\: Aligners \:\: BWA \:\: MinQuality' above.

    **Options \:\: Aligners \:\: Bowtie2 \:\: FilterUnmappedReads**
        .. literalinclude:: makefile.yaml
            :language: yaml
            :linenos:
            :lineno-start: 55
            :lines: 55-56

        See 'Options \:\: Aligners \:\: BWA \:\: FilterUnmappedReads' above.

    **Options \:\: Aligners \:\: BWA \:\: \***

        Additional command-line options may be specified for Bowtie2, as described in the "Specifying command-line options" section above. Please refer to the `Bowtie2 documentation`_ for more information about available command-line options.


mapDamage options
^^^^^^^^^^^^^^^^^
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 67
        :lines: 67-71

    This subsection is used to specify options for mapDamage2.0, when plotting *post-mortem* DNA damage, when building models of the *post-mortem* damage, and when rescaling quality scores to account for this damage. In order to enable plotting, modeling, or rescaling of quality scores, please see the 'mapDamage' option in the 'Features' section below.

    .. note::
        It may be worthwhile to tweak mapDamage parameters before building a model of *post-mortem* DNA damage; this may be accomplished by running the pipeline without rescaling, running with the 'mapDamage' feature set to 'plot' (with or without quotes), inspecting the plots generated per-library, and then tweaking parameters as appropriate, before setting 'mapDamage' to 'model' (with or without quotes).

        Should you wish to change the modeling and rescaling parameters, after having already run the pipeline with rescaling enabled, simply remove the mapDamage files generated for the relevant libraries (see the :ref:`bam_filestructure` section).

    .. warning::
        Rescaling requires a certain minimum number of C>T and G>A substitutions, before it is possible to construct a model of *post-mortem* DNA damage. If mapDamage fails with an error indicating that "DNA damage levels are too low", then it is necessary to disable rescaling for that library to continue.


**Options \:\: mapDamage :: --downsample**
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 69
        :lines: 69-71

    By default the BAM pipeline only samples 100k reads for use in constructing mapDamage plots; in our experience, this is sufficient for accurate plots and models. If no downsampling is to be done, this value can set to 0 to disable this features::

        --downsample: 100000   # Sample 100 thousand reads
        --downsample: 1000000  # Sample 1 million reads
        --downsample: 0        # No downsampling


**Options \:\: mapDamage :: \***

    Additional command-line options may be supplied to mapDamage, just like the `--downsample` parameter, as described in the "Specifying command-line options" section above. These are used during plotting and rescaling (if enabled).


Excluding read-types
^^^^^^^^^^^^^^^^^^^^
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 73
        :lines: 73-87

    During the adapter-trimming and read-merging step, AdapterRemoval will generate a selection of different read types. This option allows certain read-types to be excluded from further analyses. In particular, it may be useful to exclude non-collapsed (paired and singleton) reads when processing (ancient) DNA for which only short inserts is expected, since this may help exclude exogenous DNA. The following read types are currently recognized:

    *Single*
        Single-end reads; these are the (trimmed) reads generated from supplying single-end FASTQ files to the pipeline.

    *Paired*
        Paired-end reads; these are the (trimmed) reads generated from supplying paired-end FASTQ files to the pipeline, but covering only the subset of paired reads for which *both* mates were retained, and which were not merged into a single read (if --collapse is set for AdapterRemoval).

    *Singleton*
        Paired-end reads; these are (trimmed) reads generated from supplying paired-end FASTQ files to the pipeline, but covering only those reads in which one of the two mates were discarded due to either the `--maxns`, the `--minlength`, or the `--maxlength` options supplied to AdapterRemoval. Consequently, these reads are mapped and PCR-duplicate filtered in single-end mode.

    *Collapsed*
        Paired-end reads, for which the sequences overlapped, and which were consequently merged by AdapterRemoval into a single sequence (enabled by the --collapse command-line option). These sequences are expected to represent the complete insert, and while they are mapped in single-end mode, PCR duplicate filtering is carried out in a manner that treats these as paired reads. Note that all collapsed reads are tagged by prefixing the read name with 'M\_'.

    *CollapsedTruncated*
        Paired-end reads (like *Collapsed*), which were trimmed due to the `--trimqualities` or the `--trimns` command-line options supplied to AdapterRemoval. Consequently, and as these sequences represent the entire insert, these reads are mapped and PCR-duplicate filtered in single-end mode. Note that all collapsed, truncated reads are tagged by prefixing the read name with 'MT\_'.

    To enable / disable exclusion of a read type, set the value for the appropriate type to 'yes' or 'no' (without quotes)::

        Singleton: no   # Singleton reads are NOT excluded
        Singleton: yes  # Singleton reads are excluded


Optional features
^^^^^^^^^^^^^^^^^
    .. literalinclude:: makefile.yaml
        :language: yaml
        :linenos:
        :lineno-start: 104
        :lines: 89-107

    This section lists several optional features, in particular determining which BAM files and which summary statistics are generated when running the pipeline. Currently, the following options are available:


    *PCRDuplicates*
        This option determines how the BAM pipeline handles PCR duplicates following the mapping of trimmed reads. At present, 3 possible options are available. The first option is 'filter', which corresponds to running Picard MarkDuplicates and 'paleomix rmdup_collapsed' on the input files, and removing any read determined to be a PCR duplicate; the second option 'mark' functions like the 'filter' option, except that reads are not removed from the output, but instead the read flag is marked using the 0x400 bit (see the `SAM/BAM specification`_ for more information), in order to allow down-stream tools to identify these as duplicates. The final option is 'no' (without quotes), in which case no PCR duplicate detection / filtering is carried out on the aligned reads, useful for data generated using amplification free sequencing.


    *mapDamage*
        The 'mapDamage' option accepts four possible values: 'no', 'plot', 'model', and 'rescale'. By default value ('plot'), will cause mapDamage to be run in order to generate simple plots of the *post-mortem* DNA damage rates, as well as base composition plots, and more. If set to 'model', mapDamage will firstly generate the plots described for 'plot', but also construct models of DNA damage parameters, as described in [Jonsson2013]_. Note that a minimum amount of DNA damage is required to be present in order to build these models. If the option is set to 'rescale', both plots and models will be constructed using mapDamage, and in addition, the quality scores of bases will be down-graded based on how likely they are to represent *post-mortem* DNA damage (see above).

    *Coverage*
        If enabled, a table summarizing the number of hits, the number of aligned bases, bases inserted, and bases deleted, as well as the mean coverage, is generated for each reference sequence, stratified by sample, library, and contig.

    *Depths*
        If enabled, a table containing a histogram of the depth of coverage, ranging from 0 to 200, is generated for each reference sequence, stratified by sample, library, and contig. These files may further be used by the Phylogenetic pipeline, in order to automatically select a maximum read depth during SNP calling (see the :ref:`phylo_usage` section for more information).

    *Summary*
        If enabled, a single summary table will be generated per target, containing information about the number of reads processed, hits and fraction of PCR duplicates (per prefix and per library), and much more.

    For a description of where files are placed, refer to the :ref:`bam_filestructure` section. It is possible to run the BAM pipeline without any of these options enabled, and this may be useful in certain cases (if only the statistics or per-library BAMs are needed). To enable / disable a features, set the value for that feature to 'yes' or 'no' (without quotes)::

        Summary: no   # Do NOT generate a per-target summary table
        Summary: yes  # Generate a per-target summary table


Prefixes
--------

.. literalinclude:: makefile.yaml
    :language: yaml
    :linenos:
    :lineno-start: 110
    :lines: 110-126


Reference genomes used for mapping are specified by listing these (one or more) in the 'Prefixes' section. Each reference genome is associated with a name (used in summary statistics and as part of the resulting filenames), and the path to a FASTA file which contains the reference genome. Several other options are also available, but only the name and the 'Path' value are required, as shown here for several examples::

    #  Map of prefixes by name, each having a Path key, which specifies the
    # location of the BWA/Bowtie2 index, and optional label, and an option
    # set of regions for which additional statistics are produced.
    Prefixes:
      # Name of the prefix; is used as part of the output filenames
      MyPrefix1:
        # Path to FASTA file containing reference genome; must end with '.fasta'
        Path: /path/to/genomes/file_1.fasta
      MyPrefix2:
        Path: /path/to/genomes/file_2.fasta
      MyPrefix3:
        Path: /path/to/genomes/AE008922_1.fasta

Each sample in the makefile is mapped against each prefix, and BAM files are generated according to the enabled 'Features' (see above). In addition to the path, it is possible to specify `RegionsOfInterest`, which are described below.


Regions of interest
^^^^^^^^^^^^^^^^^^^

It is possible to specify one or more "regions of interest" for a particular reference genome. Doing so results in the production of coverage and depth tables being generated for those regions (if these features are enabled, see above), as well as additional information in the summary table (if enabled).

Such regions are specified using a BED file containing one or more regions; in particular, the first three columns (name, 0-based start coordinate, and 1-based end coordinate) are required, with the 4th column (the name) being optional. Strand information (the 6th column) is not used, but must still be valid according to the BED format.

Statistics are merged by the names specified in the BED file, or by the contig names if no names were specified. Thus, it is important to insure that names are unique if individual statistics are desired for every region.

Specifying regions of interest is accomplished by providing a name and a path for each set of regions of interest under the `RegionOfInterest` section for a given prefix::

    # Produce additional coverage / depth statistics for a set of
    # regions defined in a BED file; if no names are specified for the
    # BED records, results are named after the chromosome / contig.
    RegionsOfInterest:
      MyRegions: /path/to/my_regions.bed
      MyOtherRegions: /path/to/my_other_regions.bed

The following is a simple example of such a BED file, for an alignment against the rCRS (`NC_012920.1`_)::

    NC_012920_1 3306    4262    region_a
    NC_012920_1 4469    5510    region_b
    NC_012920_1 5903    7442    region_a

In this case, the resulting tables will contain information about two different regions, namely `region_a` (2495 bp, resulting from merging the two individual regions specified), and `region_b` (1041 bp). The order of lines in this file does not matter.


Adding multiple prefixes
^^^^^^^^^^^^^^^^^^^^^^^^

In cases where it is necessary to map samples against a large number of reference genomes, it may become impractical to add these to the makefile by hand. It is therefore possible to specify the location of the reference genomes via a path containing wild-cards, and letting the BAM pipeline collect these automatically. For the following example, we assume that we have a folder at '/path/to/genomes', which contains our reference genomes:

.. code-block:: bash

    $ ls /path/to/genomes
    AE000516_2.fasta
    AE004091_2.fasta
    AE008922_1.fasta
    AE008923_1.fasta

To automatically add these four reference genomes to the makefile, we would add a prefix as follows::

    #  Map of prefixes by name, each having a Path key, which specifies the
    # location of the BWA/Bowtie2 index, and optional label, and an option
    # set of regions for which additional statistics are produced.
    Prefixes:
      # Name of the prefix; is used as part of the output filenames
      MyGenomes*:
        # Path to .fasta file containing a set of reference sequences.
        Path: /path/to/genomes/*.fasta

There are two components to this, namely the name of the pseudo-prefix which *must* end with a star (\*), and the path which may contain one or more wild-cards. If the prefix name does not end with a star, the BAM pipeline will simply treat the path as a regular path. In this particular case, the BAM pipeline will perform the equivalent of 'ls /path/to/genomes/\*.fasta', and then add each file it has located using the filename without extensions as the name of the prefix. In other words, the above is equivalent to the following::

    #  Map of prefixes by name, each having a Path key, which specifies the
    # location of the BWA/Bowtie2 index, and optional label, and an option
    # set of regions for which additional statistics are produced.
    Prefixes:
      # Name of the prefix; is used as part of the output filenames
      AE000516_2:
        Path: /path/to/genomes/AE000516_2.fasta
      AE004091_2:
        Path: /path/to/genomes/AE004091_2.fasta
      AE008922_1:
        Path: /path/to/genomes/AE008922_1.fasta
      AE008923_1:
        Path: /path/to/genomes/AE008923_1.fasta

.. note::
    The name provided for the pseudo-prefix (here 'MyGenomes') is not used by the pipeline, and can instead be used to document the nature of the files being included.

.. warning::
    Just like regular prefixes, it is required that the filename of the reference genome ends with '.fasta'. However, the pipeline will attempt to add *any* file found using the provided path with wildcards, and care should therefore be taken to avoid including non-FASTA files. For example, if the path '/path/to/genomes/\*' was used instead of '/path/to/genomes/\*.fasta', this would cause the pipeline to abort due to the inclusion of (for example) non-FASTA index files generated at this location by the pipeline itself.


Targets
-------
.. literalinclude:: makefile.yaml
    :language: yaml
    :linenos:
    :lineno-start: 129
    :lines: 129-

In the BAM pipeline, the term 'Target' is used to refer not to a particular sample (though in typical usage a target includes just one sample), but rather one or more samples to processed together to generate a BAM file per prefix (see above). A sample included in a target may likewise contain one or more libraries, for each of which one or more sets of FASTQ reads are specified.

The following simplified example, derived from the makefile constructed as part of :ref:`bam_usage` section exemplifies this:

.. code-block:: yaml
    :linenos:

    # Target name; all output files uses this name as a prefix
    MyFilename:
      # Sample name; used to tag data for segregation in downstream analyses
      MySample:
        # Library name; used to tag data for segregation in downstream analyses
        TGCTCA:
          # Lane / run names and paths to FASTQ files
          Lane_1: data/TGCTCA_L1_*.fastq.gz
          Lane_2: data/TGCTCA_L2_R{Pair}_*.fastq.gz


*Target name*
    The first top section of this target (line 2, 'MyFilename') constitute the target name. This name is used as part of summary statistics and, more importantly, determined the first part of name of files generated as part of the processing of data specified for this target. Thus, in this example all files and folders generated during the processing of this target will start with 'MyFilename'; for example, the summary table normally generated from running the pipeline will be placed in a file named 'MyFilename.summary'.

*Sample name*
    The subsections listed in the 'Target' section (line 4, 'MySample') constitute the (biological) samples included in this target; in the vast majority of analyses, you will have only a single sample per target, and in that case it is considered good practice to use the same name for both the target and the sample. A single target can, however, contain any number of samples, the data for which are tagged according to the names given in the makefile, using the SAM/BAM readgroup ('RG') tags.

*Library name*
    The subsections listed in the 'Sample' section (line 6, 'TGCTCA') constitute the sequencing libraries constructed during the extraction and library building for the current sample. For modern samples, there is typically only a single library per sample, but more complex sequencing projects (modern and ancient) may involve any number of libraries constructed from one or more extracts. It is very important that libraries be listed correctly (see below).

    .. warning::
        Note that the BAM pipeline imposes the restriction that each library name specified for a target must be unique, even if these are located in to different samples. This restriction may be removed in future versions of PALEOMIX.

*Lane name*
    The subsections of each library are used to specify the names of individual lanes. Each sample may include any number of lanes and each lane may include any number of FASTQ files. Paths listed here may include wildcards (\*) or the special value `{Pair}`, which is used to indicate that the lane is paired end. During runtime, `{Pair}` is replaced with `1` and `2`, before searching using wildcards, and the pipeline expects the resulting filenames to correspond to mate 1 and mate 2 reads respectively.

In addition to these target (sub)sections, it is possible to specify 'Options' for individual targets, samples, and libraries, similarly to how this is done globally at the top of the makefile. This is described below.

.. warning::
    It is very important that lanes are assigned to their corresponding libraries in the makefile; while it is possible to simply record every sequencing run / lane under a single library and run the pipeline like that, this will result in several unintended side effects: Firstly, the BAM pipeline uses the library information to ensure that PCR duplicates are filtered correctly. Wrongly grouping together lanes will result either in the loss of sequences which are not, in fact, PCR duplicates, while wrongly splitting a library into multiple entries will result in PCR duplicates not being correctly identified across these. Furthermore, mapDamage analyses make use of this information to carry out various analyses on a per-library basis, which may similarly be negatively impacted by incorrect specification of libraries.


Including already trimmed reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some cases it is useful to include FASTQ reads that has already been trimmed for adapter sequences. While this is not recommended in general, as it may introduce systematic bias if some data has been processed differently than the remaining FASTQ reads, the BAM pipeline makes it simple to incorporate both 'raw' and trimmed FASTQ reads, and to ensure that these integrate in the pipeline.

To include already trimmed reads, these are specified as values belonging to a lane, using the same names for read-types as in the 'ExcludeReads' option (see above). The following minimal example demonstrates this:

.. code-block:: yaml
    :linenos:

    MyFilename:
      MySample:
        ACGATA:
          # Regular lane, containing reads that are not already trimmed
          Lane_1: data/ACGATA_L1_R{Pair}_*.fastq.gz

          # Lane containing pre-trimmed reads of each type
          Lane_2:
            # Single-end reads
            Single: /path/to/single_end_reads.fastq.gz

            # Paired-end reads where one mate has been discarded
            Singleton: /path/to/singleton_reads.fastq.gz

            # Paired end reads; note that the {Pair} key is required,
            # just like with raw, paired-end reads
            Paired: /path/to/paired_end_{Pair}.fastq.gz

            # Paired-end reads merged into a single sequence
            Collapsed: /path/to/collapsed.fastq.gz

            # Paired-end reads merged into a single sequence, and then truncated
            CollapsedTruncated: /path/to/collapsed_truncated.fastq.gz

The above examples show how each type of reads are to be listed, but it is not necessary to specify more than a single type of pre-trimmed reads in the makefile.

.. note::
    Including already trimmed reads currently result in the absence of some summary statistics in the .summary file, namely the number of raw reads, as well as trimming statistics, since the BAM pipeline currently relies on AdapterRemoval to collect these statistics.

Overriding global settings
^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the main `Options` section, it is possible to override options at a Target, Sample, and Library level. This allows, for example, that different adapter sequences be specified for each library generated for a sample, or options that should only be applied to a particular sample among several included in a makefile. The following demonstration uses the makefile constructed as part of :ref:`bam_usage` section as the base:

.. code-block:: yaml
    :linenos:
    :emphasize-lines: 2-7, 10-14, 20-23

    MyFilename:
      # These options apply to all samples with this filename
      Options:
        # In this example, we override the default adapter sequences
        AdapterRemoval:
          --adapter1: AGATCGGAAGAGC
          --adapter2: AGATCGGAAGAGC

      MySample:
        # These options apply to libraries 'ACGATA', 'GCTCTG', and 'TGCTCA'
        Options:
           # In this example, we assume that FASTQ files for our libraries
           # include Phred quality scores encoded with offset 64.
           QualityOffset: 64

        ACGATA:
          Lane_1: data/ACGATA_L1_R{Pair}_*.fastq.gz

        GCTCTG:
          # These options apply to 'Lane_1' in the 'GCTCTG' library
          Options:
            # It is possible to override options we have previously overridden
            QualityOffset: 33

          Lane_1: data/GCTCTG_L1_*.fastq.gz

        TGCTCA:
          Lane_1: data/TGCTCA_L1_*.fastq.gz
          Lane_2: data/TGCTCA_L2_R{Pair}_*.fastq.gz


In this example, we have overwritten options at 3 places:

* The first place (lines 2 - 7) will be applied to *all* samples, libraries, and lanes in this target, unless subsequently overridden. In this example, we have set a new pair of adapter sequences, which we wish to use for these data.

* The second place (lines 10 - 14) are applied to the sample 'MySample' that we have included in this target, and consequently applies to all libraries specified for this sample ('ACGATA', 'GCTCTG', and 'TGCTCA'). In most cases you will only have a single sample, and so it will not make a difference whether or not you override options for the entire target (e.g. lines 3 - 8), or just for that sample (e.g. lines 11-15).

* Finally, the third place (lines 20 - 23) demonstrate how options can be overridden for a particular library. In this example, we have chosen to override an option (for this library only!) we previously overrode for that sample (the 'QualityOffset' option).

.. note:: It currently not possible to override options for a single lane, it is only possible to override options for all lanes in a library.

.. warning::
    Only the `mapDamage` and the `PCRDuplicates` features can be overridden for individual targets, samples, or libraries. Over features can be overridden per target, but not per sample or per library.


.. _AdapterRemoval documentation: https://github.com/MikkelSchubert/adapterremoval
.. _Bowtie2 documentation: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _NC_012920.1: http://www.ncbi.nlm.nih.gov/nuccore/251831106
.. _Phred quality-scores: https://en.wikipedia.org/wiki/FASTQ_format#Quality
.. _SAM/BAM specification: http://samtools.sourceforge.net/SAM1.pdf
.. _seqtk: https://github.com/lh3/seqtk
