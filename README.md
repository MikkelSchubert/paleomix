Installation instructions for BAM pipeline
==========================================
The following instructions will use ~/install as the installation directory for the pipeline and required tools, but any directory may be used.

The pipeline currently requires Python v2.6 or 2.7, in addition to a few modules, as well as [git](http://git-scm.com/) and a few other tools (see below).

1. Installing pipeline
----------------------
To install a copy of the pipeline, create a git clone as follows:

    $ mkdir -p ~/install
    $ cd ~/install
    $ git clone https://github.com/MikkelSchubert/pypeline.git
    <login using your GitHub account information>

The pypeline can be installed in a couple of ways:

The least invasive installation is accomplished by updating your PYTHONPATH and PATH variables as follows (for bash users):

    $ echo "export PYTHONPATH=\$PYTHONPATH:~/install/pypeline" >> ~/.bashrc
    $ echo "export PATH=\$PATH:~/install/pypeline/bin" >> ~/.bashrc

Alternatively, you can make a local install using the included 'setup.py' script:

    $ python setup.py install --user

Finally, the pipeline can be installed for all users using the same script:

    $ sudo python setup.py install

Regardless of installation method, you may need to re-login for these changes to take effect.


2. Installing required modules
------------------------------
The pipeline requires [Cython 0.18+](http://www.cython.org/), [Pysam v0.7.4+](https://code.google.com/p/pysam/), and [PyYAML v3.1+](http://www.pyyaml.org/). Please note that Cython MUST be installed before Pysam. Each package can be installed using the following command, in the root of the source folders:

    $ python setup.py install --user

For example

    $ wget http://pyyaml.org/download/pyyaml/PyYAML-3.10.tar.gz
    $ tar xvzf PyYAML-3.10.tar.gz
    $ cd PyYAML-3.10
    $ python setup.py install --user


2. Installing required applications
-----------------------------------
The BAM pipeline requires [SAMTools](http://samtools.sourceforge.net) v0.1.18+, [AdapterRemoval](http://code.google.com/p/adapterremoval/) v1.4+ (1.5+ recommended for new projects, and please cite [Lindgreen 2013](http://www.ncbi.nlm.nih.gov/pubmed/22748135)!), and [BWA](http://bio-bwa.sourceforge.net/) v0.5.9+ / v0.7.5+ or [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 2.0.0+.  Only the aligner (BWA/Bowtie2) that is to be used (selected in the individual makefiles) needs to be installed.

In addition, the pipeline requires the [Picard](http://picard.sourceforge.net/) v1.82+ jars. These instructions assume that these jars have been placed in ~/install/jar_root, but any location may be used. This folder is specified using the --jar-root parameter (see below).

If the user wishes to perform local realignment around indels (recommended), then [GATK](http://www.broadinstitute.org/gatk/) is required. The GATK jar file MUST be placed in the same folder as the Picard jars. Note that if alignments are to be carried out against the human nuclear genome, chromosome MUST be ordered by their number for GATK to work! See the GATK website for more information / files.


3. Setting default options
--------------------------
Every command-line option for the pipeline, which has a default value, can be specified using the ~/.pypeline.conf file, by removing the first two dashes, and replacing other dashes with underscores. It is recommended to create a ~/.pypeline.conf file, containing at least values for --jar-root and for --temp-root. For example:

    $ cat ~/.pypeline.conf
    [Defaults]
    temp_root = ~/temp/bam_pipeline
    jar_root  = ~/install/jar_root

The jar\_root (--jar-root) and temp\_root (--temp-root) parameters default to the value shown above, but other paths may be used as desired.

WARNING: Using /tmp for --temp-root is NOT recommended, since the pipeline generates many large files. If care is not taken, it is very easy to fill up /tmp, causing significant problems for EVERYONE using the server.

It is recommended to use a folder on the same partition as the destination folder, to reduce disk usage, as all generated files will have to be moved from the --temp\_root folder to the --destination folder during the operation of the pipeline. It is furthermore recommended to use a separate subfolder (ie. the "bam_pipeline" folder in the example above), since failed runs may leave many temporary files behind in order to facilitate debugging.


Creating a makefile
===================
The 'bam_pipeline mkfile' command can be used to create a makefile template. In addition to the various options for the programs run during the pipeline (documented in the template), it is nessesary to specify the name and location of target sequences (termed Prefixes) and each sample to be processed (grouped under one or more targets). An example project, based on synthetic reads genereted from the rCRS, is included in the distribution under examples/bam_pipeline.

If 'SampleSheet.csv' files are present, these may be specified as arguments to 'bam_pipeline mkfile', at which point the information found in the file is incorported into the makefile template. The quality of the resulting makefile will depend heavily on the quality of the information in the 'SampleSheet.csv' files, so it is recommended that the output be inspected carefully:

    $ bam_pipeline mkfile reads/SampleSheet.csv > Project.makefile


1. Prefixes
-----------
Prefixes are built using either the BWA command "bwa index /path/to/sequence.fasta" or the Bowtie2 command "bowtie2-build /path/to/sequence.fasta /path/to/sequence.fasta" (the path is intentionally repeated). For example, assuming that we wish to align against a FASTA file located at 'prefixes/genome.fa', we may index the file as follows:

    $ bwa index prefixes/genome.fa

Adding this to the makefile is accomplished by adding the following files:

    Prefixes:
      my_genome:
        Path: prefixes/genome.fa

The name of the prefix (here 'my_genome') will be used to name the resulting files and in resulting tables that are generated. Typical names include 'hg19', 'EquCab20', and other standard abbrivations for reference genomes, accession numbers, and the like. Multiple prefixes can be specified, but each name MUST be unique:

    Prefixes:
      my_genome:
        Path: prefixes/genome.fa
      another_genome:
        Path: prefixes/genome2.fa


Additional options for the prefixes are documented below.


2. Target information
---------------------
Each makefile must contain one or more targets, which are described by 4 properties, and a path:

1. A unique target name, used as a prefix for the resulting output files.
2. One or more samples for a target, used to set the SM readgroup tag in BAM files. For simple projects with only one sample, this can simply be the same name as the target.
3. One or more libraries for a sample, used to set the LB readgroup tag in BAM files. It is recommended to use a name that allows for the subsequent identification of the physical library used to generate these reads, for exampel the index used during multiple plexing. Note that duplicates are removed per library!
4. One or more lanes for a library, used to set the PU readgroup tag in BAM files.

In addition to this information, a path to the FASTQ files are required. This path may contain wild-cards, if multiple files were generated for a given lane.

For example, assuming that we have a sample which we have given the barcode 'ED209', from which we have run one lane (#2) of the library using the index 'ACGTTA', which we located in 'reads/ED209_ACGTTA_L002_R1_*.fastq.gz' we might use the following structure:

    ED209:
      ED209:
	    ACGTTA:
		   Lane01: 'reads/ED209_ACGTTA_L002_R1_*.fastq.gz'

Paired-ended data is expected to be split into two sets of files, one for each mate, and are specified by adding the '{Pair}' value to the path, which is replaced by '1' and by '2' during setup. For example, mates are typically specified using 'R1' and 'R2' for Illumina data, so the example above for PE data would become:

    ED209:
      ED209:
	    ACGTTA:
		   Lane01: 'reads/ED209_ACGTTA_L002_R{Pair}_*.fastq.gz'

As noted above, any number of targets, samples, libraries and lanes may be specified:

    Target_1:
	  Sample_1:
	    Library_1:
		  Lane_1: ...
		  Lane_2: ...
	    Library_2:
		  ...
	  Sample_2:
	    ...
	Target_2:
	  ...

3. Additional options for prefixes
----------------------------------
Each prefix may be assigned a label (either 'mitochondrial' or 'nuclear'), which is useful if these two are aligned sepearately. If both labels are specified, the summary file that is generated will contain information aggregating the information for each of these prefixes under the name 'endogenous'. Each label may be specified any number of times. For example:

    Prefixes:
      hg19:
        Path:  prefixes/hg19.fa
	Label: nuclear
      chrM:
        Path:  prefixes/rCRS.fa
	Label: mitochondrial

In addition, it is possible to specify multiple prefixes at once using wildcards. To enable this option, postfix the name of the prefix with a wildcard ('*') and add wildcards to the Path. This is helpful if sequences are to be aligned against many genomes. For example, if the folder 'indices' contains 'hg19.fa' and 'rCRS.fa', the following settings would be equivalent to the above, except that no Labels are set:

    Prefixes:
      human*:
        Path: prefixes/*.fa

If a Label is specified, it is applied to all prefixes that are found.

Finally, areas of interest may be specified for which coverage and depth information is to be calculated. Each area has a name (used in the output files), and points to a bedfile containing the relevant regions of the genome. If no names are given in the BED file, the intervals are merged by contig, and named after the contig with a wildcard ("*") appended. The following example demonstrates how this may be acomplished, using a bedfile located at 'prefixes/hg19.exome.bed':

    Prefixes:
      hg19:
        Path:  prefixes/hg19.fa
        AreasOfInterest:
          Exome: prefixes/hg19.exome.bed


4. Additional options for samples
----------------------------------
By default, all files specified for lanes are trimmed using AdapterRemoval. If reads have already been trimmed, these can be specified as follows, using the same read-types as in the 'ExcludeReads' options. For example:

    ED209:
      ED209:
	    ACGTTA:
		   Lane01:
		    Single: reads/trimmed.fa.gz # Path to SE reads
			Paired: reads/trimmed.{Pair}.fa.gz # Path to PE reads (must have {Pair} component)

If reads have already been mapped, the BAM can be incorported into the project, in which case the BAM is retagged using the specified properties, and merged into the final BAM file. This is done by specifying the name of the prefix instead of a read-type. These are assumed to contain SE and/or PE reads, and not collapsed reads. For example, assuming that we are using the prefix 'hg19':

    ED209:
      ED209:
	    ACGTTA:
		   Lane01:
		    hg19: bams/old.bam




Citations
========================

 * Langmead B, Salzberg SL. "Fast gapped-read alignment with Bowtie 2" Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.
 * Li H, Durbin R. "Fast and accurate short read alignment with Burrows-Wheeler transform" Bioinformatics. 2009 Jul 15;25(14):1754-60. doi:10.1093/bioinformatics/btp324.
 * Li H, *et al*. "The Sequence Alignment/Map format and SAMtools" Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352.
 * Lindgreen S. "AdapterRemoval: easy cleaning of next-generation sequencing reads" BMC Res Notes. 2012 Jul 2;5:337. doi: 10.1186/1756-0500-5-337.
 * McKenna A, *et al*. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data" Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110.
