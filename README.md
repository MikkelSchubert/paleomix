BAM / Phylogenetic pipelines
==========================================

TODO



---------------------
# 1. Installation instructions
The following instructions will use ~/install as the installation
directory for the pipelines and required tools, but any directory may
be used. The two pipelines (BAM and phylogenetic) have their own set
of requirements, as well as a common set of requirements.

Note that while it is possible to install and run the pipelines on
MacOSX (e.g. Mountain Lion), this is not recommended. Users who wish
to do so are strongly encouraged to use the
[HomeBrew](http://brew.sh/) package manager.


## 1.1. Common instructions
The following instructions apply to both pipelines included in the
repository. The pipelines require [Python 2.7](http://python.org), in
addition to a few modules. [Git](http://git-scm.com/) is recommended
to ease upgrading of the pipeline (see below).


### 1.1.1. Obtaining the pypeline source-code
To install a copy of the pipeline, create a git clone as follows:

    $ mkdir -p ~/install
	$ cd ~/install
	$ git clone https://github.com/MikkelSchubert/pypeline.git

or download a ZIP archive containing the repository:

    $ mkdir -p ~/install
	$ cd ~/install
	$ wget	https://github.com/MikkelSchubert/pypeline/archive/master.zip
	$ mv pypeline-master pypeline


### 1.1.2 Installing the pipeline
The least invasive installation is accomplished by updating your
PYTHONPATH and PATH variables as follows (for Bash users):

    $ echo "export PYTHONPATH=\$PYTHONPATH:~/install/pypeline" >> ~/.bashrc
    $ echo "export PATH=\$PATH:~/install/pypeline/bin" >> ~/.bashrc

Alternatively, you can make a local install using the included
'setup.py' script:

    $ cd ~/install/pypeline
	$ python setup.py install --user

Finally, the pipeline can be installed for all users using the same
script:

    $ cd ~/install/pypeline
	$ sudo python setup.py install

Regardless of installation method, you may need to start a new Bash
session for these changes to take effect.


### 1.1.3. Installing required modules The pipeline requires
[Pysam](https://code.google.com/p/pysam/)  v0.7.5+:

    $ curl -O https://pysam.googlecode.com/files/pysam-0.7.5.tar.gz
	$ tar xvzf pysam-0.7.5.tar.gz
	$ cd pysam-0.7.5
	$ python setup.py install --user


## 1.2 Installing programs required by the BAM pipeline

The BAM pipeline requires the following programs:

 * [AdapterRemoval](http://code.google.com/p/adapterremoval/) v1.5+
   [[Lindgreen 2013](http://www.ncbi.nlm.nih.gov/pubmed/22748135)]
 * [BEDTools](https://code.google.com/p/bedtools/) v2.16.0+
   [[Quinlan and Hall 2010](http://www.ncbi.nlm.nih.gov/pubmed/20110278)]
 * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) 2.1.0+
   [[Langmead and Salzberg 2012](http://www.ncbi.nlm.nih.gov/pubmed/22388286)]
 * [BWA](http://bio-bwa.sourceforge.net/) v0.5.9+, or 0.6.2+, or
   v0.7.5+ [[Li and Durbin 2009](http://www.ncbi.nlm.nih.gov/pubmed/19451168)]
 * [SAMTools](http://samtools.sourceforge.net) v0.1.18+
   [[Li et al 2009](http://www.ncbi.nlm.nih.gov/pubmed/19505943)]

Either BWA and Bowtie2 needs to be installed, corresponding to the
program that the user enables in the makefiles (see below).

  * [mapDamage](http://ginolhac.github.io/mapDamage/) v2.0.2+
    [[Hákon et al 2013](http://www.ncbi.nlm.nih.gov/pubmed/23613487)]

Note that the GNU Scientific Library (GSL) and the R packages listed
in the mapDamage instructions are required for rescaling support;
inline, gam, Rcpp, RcppGSL and ggplot2 (>=0.9.2). Use the following
commands to verify that GSL / R packages have been correctly
installed:

    $ mapDamage --check-R-packages
	$ gsl-config

In addition, the JARs for the following programs are required:
 * [Picard tools](http://picard.sourceforge.net/) v1.82+
 * [GATK](http://www.broadinstitute.org/gatk/) (any version)
   [[McKenna et al 2010](http://www.ncbi.nlm.nih.gov/pubmed/20644199)]


The Picard Tools JARs (.jar) files are expected to be located in
~/install/jar_root/ by default, but this behaviour may be changed
using either the --jar-root command-line option, or the global
configuration file (see below).

The GATK JAR is only required if the user wishes to carry out local
realignment near indels (recommended), and is expected to be placed in
the same folder as the Picard Tools JARs.

**WARNING:** If alignments are to be carried out against the human
nuclear genome, chromosomes MUST be ordered by their number for GATK
to work!  See the GATK
[FAQ](http://www.broadinstitute.org/gatk/guide/topic?name=faqs#id1204)
for more information.


## 1.3 Installing programs required by the phylogenetic pipeline

The following programs are required for various parts of the
pipeline. The version numbers given below specify the minimum version
required to run the pipeline.

**Genotyping:**

 * [BEDTools](https://code.google.com/p/bedtools/) v2.16.0+
   [[Quinlan and Hall 2010](http://www.ncbi.nlm.nih.gov/pubmed/20110278)]
 * [SAMTools](http://samtools.sourceforge.net) v0.1.18+
   [[Li et al 2009](http://www.ncbi.nlm.nih.gov/pubmed/19505943)]
 * [Tabix](http://samtools.sourceforge.net/) v0.2.5

Both the 'tabix' and the 'bgzip' executable from the Tabix package
must be installed.

**Multiple Sequence Alignment:**

* [MAFFT v7+](http://mafft.cbrc.jp/alignment/software/)
  [[Katoh and Standley 2013](http://www.ncbi.nlm.nih.gov/pubmed/23329690)]

Note that the pipeline requires that the algorithm-specific MAFFT
commands (e.g. 'mafft-ginsi', 'mafft-fftnsi'). These are automatically
created by the 'make install' command.

**Phylogenetic Inference:**

* [RAxML v7.3.2+](https://github.com/stamatak/standard-RAxML)
  [[Stamatakis 2006](http://www.ncbi.nlm.nih.gov/pubmed/16928733)]
* [ExaML v1.0.5+, parser v1.0.2+](https://github.com/stamatak/ExaML)

The pipeline expects a single-threaded binary named 'raxmlHPC' for
RAxML. The pipeline expects the ExaML binary to be named 'examl', and
the parser binary to be named 'examlParser'. Compiling and running
ExaML requires an MPI implementation (e.g. Open-MPI;
http://www.open-mpi.org/), even if ExaML is run single-threaded.

Both programs offer a variety of makefiles suited for different
server-architectures and use-cases. If in doubt, use the
Makefile.SSE3.gcc makefiles, which are compatible with most modern
systems:

    $ make -f Makefile.SSE3.gcc


## 1.4. Setting default options

Default locations of temporary files, JAR files, the maximum number of
active processes, and other options may be set globally or per host
(if an installation is shared across multiple servers). This is
accomplished by using the --write-config-file switch, which writes the
current (globally settable) values to ~/.pypeline/bam_pipeline.ini and
~/.pypeline/phylo_pipeline.ini:

    $ bam_pipeline run --write-config-file
	$ phylo_pipeline run --write-config-file

It is possible to specify which values to write the config file, by
passing them together with --write-config-file. The configuration
files contain a single “[Defaults]” section, which is applied to all
hosts on which the pipeline is run:

    [Defaults]
    max_threads = 4
    ...

It is furthermore possible to over-ride settings on a per-host basis
(e.g. to set the number of threads or the temporary directory), by
adding a section named using the host-name of the server in question,
(as determined by the “hostname” command):

    [Defaults]
    max_threads = 4
    ...
    [MyHostName]
    max_threads = 16
    ...

**WARNING**: By default, the pipeline will use
/tmp/**username**/**pipeline_name** for temporary files and
folders. As this includes all results generated by the pipeline (which
are only moved to the destination upon completion of a given command),
the value *must* be changed if /tmp folder does not contain sufficient
free space, as filling up /tmp will cause failure for all users of the
server.


## 1.5 Testing the pipelines
Example projects are included for both the BAM and the phylogenetic
pipeline. Running these projects is recommended, in order to verify that
the pipeline and required applications have been correctly installed:

*BAM pipeline:*

    $ cd ~/install/pypeline/examples/bam_pipeline/
	$ bam_pipeline run 000_makefile.yaml

*Phylogenetic pipeline:*

    $ cd ~/install/pypeline/examples/bam_pipeline/
	$ ./setup.sh
	$ phylo_pipeline genotyping+msa+phylogeny 000_makefile.yaml


---------------------
# 2. BAM pipeline usage

## 2.1 Creating a makefile
The 'bam\_pipeline mkfile' command can be used to create a makefile
template in [YAML](http://www.yaml.org) format. In addition to the
various options for the programs run during the pipeline (documented
in the template), it is necessary to specify the name and location of
target sequences (termed Prefixes) and each sample to be processed
(grouped under one or more targets).

    $ bam_pipeline mkfile > makefile.yaml

If 'SampleSheet.csv' files are present, these may be specified as
arguments to 'bam_pipeline mkfile', at which point the information
found in the file is incorporated into the makefile template. The
quality of the resulting makefile will depend heavily on the quality
of the information in the 'SampleSheet.csv' files, so it is
recommended that the output be inspected carefully:

    $ bam_pipeline mkfile reads/SampleSheet.csv > makefile.yaml

An example project, based on synthetic reads generated from the rCRS,
is included in the distribution under examples/bam_pipeline. The
pipeline is executed using the 'run' command:

    $ bam_pipeline run [OPTIONS] <MAKEFILE>


### 2.1.1 Prefixes
Prefixes are built using either the BWA command "bwa index" or the
Bowtie2 command "bowtie2-build". The pipeline expects the prefixes to
correspond to the input filename, as would be the case if (for
example) "bwa index" was run as follows:

    $ bwa index prefixes/genome.fa

The pipeline will automatically carry out indexing of the FASTA files,
and therefore requires write-access to the folder containing the FASTA
files. If this is not possible, one may simply create a folder
containing symbolic links to the original FASTA files.

Specifying which FASTA file to align sequences is accomplished as
follows:

    Prefixes:
      my_genome:
        Path: prefixes/genome.fa

The name of the prefix (here 'my_genome') will be used to name the
resulting files and in resulting tables that are generated. Typical
names include 'hg19', 'EquCab20', and other standard abbrivations for
reference genomes, accession numbers, and the like. Multiple prefixes
can be specified, but each name MUST be unique:

    Prefixes:
      my_genome:
        Path: prefixes/genome.fa
      another_genome:
        Path: prefixes/genome2.fa


Additional options for the prefixes are documented below.


### 2.1.2 Target information
Each makefile must contain one or more targets, which are described by
4 properties, and a path:

 1. A unique target name, used as a prefix for the resulting output
    files.
 2. One or more samples for a target, used to set the SM readgroup tag
    in BAM files. For simple projects with only one sample, this can
    simply be the same name as the target.
 3. One or more libraries for a sample, used to set the LB readgroup
    tag in BAM files. It is recommended to use a name that allows for
    the subsequent identification of the physical library used to
    generate these reads, for exampel the index used during multiple
    plexing. Note that duplicates are removed per library!
 4. One or more lanes for a library, used to set the PU readgroup tag
    in BAM files.

In addition to this information, a path to the FASTQ files are
required. This path may contain wild-cards, if multiple files were
generated for a given lane.

For example, assuming that we have a sample which we have given the
barcode 'ED209', from which we have run one lane (#2) of the library
using the index 'ACGTTA', which we located in
'reads/ED209_ACGTTA_L002_R1_*.fastq.gz' we might use the following
structure:

    ED209:
      ED209:
        ACGTTA:
          Lane01: 'reads/ED209_ACGTTA_L002_R1_*.fastq.gz'

Paired-ended data is expected to be split into two sets of files, one
for each mate, and are specified by adding the '{Pair}' value to the
path, which is replaced by '1' and by '2' during setup. For example,
mates are typically specified using 'R1' and 'R2' for Illumina data,
so the example above for PE data would become:

    ED209:
      ED209:
        ACGTTA:
          Lane01: 'reads/ED209_ACGTTA_L002_R{Pair}_*.fastq.gz'

As noted above, any number of targets, samples, libraries and lanes
may be specified:

    Target_1:
      Sample_1:
        Library_1:
          Lane_1: ...
          Lane_2: ...
        Library_2: ...
      Sample_2: ...
    Target_2: ...

### 2.1.3. Additional options for prefixes
Each prefix may be assigned a label (see the template for a list of
possible values), which is useful if these two are aligned separately.
If both the 'mitochondrial' and the 'nuclear' labels are specified, the
summary file that is generated will contain information aggregating the
information for each of these prefixes under the name
'endogenous'. Each label may be specified any number of times. For
example:

    Prefixes:
      hg19:
        Path: prefixes/hg19.fa
        Label: nuclear
      chrM:
        Path: prefixes/rCRS.fa
        Label: mitochondrial

In addition, it is possible to specify multiple prefixes at once using
wildcards. To enable this option, postfix the name of the prefix with
a wildcard ('*') and add wildcards to the Path. This is helpful if
sequences are to be aligned against many genomes. For example, if the
folder 'indices' contains 'hg19.fa' and 'rCRS.fa', the following
settings would be equivalent to the above, except that no Labels are
set:

    Prefixes:
      human*:
        Path: prefixes/*.fa

If a Label is specified, it is applied to all prefixes that are found.

Finally, areas of interest may be specified for which coverage and
depth information is to be calculated. Each area has a name (used in
the output files), and points to a bedfile containing the relevant
regions of the genome. If no names are given in the BED file, the
intervals are merged by contig, and named after the contig with a
wildcard ("*") appended. The following example demonstrates how this
may be accomplished, using a BED file located at
'prefixes/hg19.exome.bed':

    Prefixes:
      hg19:
        Path: prefixes/hg19.fa
        AreasOfInterest:
          Exome: prefixes/hg19.exome.bed

### 2.1.4. Alternative input files
By default, all files specified for
lanes are trimmed using AdapterRemoval. If reads have already been
trimmed, these can be specified as follows, using the same read-types
as in the 'ExcludeReads' options. For example:

    ED209:
      ED209:
        ACGTTA:
          Lane01: Single: reads/trimmed.fa.gz # Path to SE reads
          Paired: reads/trimmed.{Pair}.fa.gz # Path to PE reads (must have {Pair} component)

If reads have already been mapped, the BAM can be incorporated into
the project, in which case the BAM is retagged using the specified
properties, and merged into the final BAM file. This is done by
specifying the name of the prefix instead of a read-type. These are
assumed to contain SE and/or PE reads, and not collapsed reads. For
example, assuming that we are using the prefix 'hg19':

    ED209:
      ED209:
        ACGTTA:
          Lane01:
            hg19: bams/old.bam


---------------------
# 3. Phylogenetic pipeline usage

## 3.1 Creating a makefile
The 'phylo\_pipeline mkfile' command can be
used to create a makefile template, as with the 'bam\_pipeline mkfile'
command (see 3.1). This makefile is used to specify the samples,
regions of interest (to be analysed), and options for the various
programs.

    $ phylo_pipeline mkfile > makefile.yaml

An example project, based on synthetic reads is included in
examples/phylo_pipeline. To run it, first execute the 'setup.sh'
script in the alignment subfolder, and run the BAM pipeline on the
makefile in that folder.

Note that filenames are not specified explicitly with this pipeline,
but are instead inferred from the names of samples, prefixes, etc. as
described below.

To execute the pipeline, a command corresponding to the step to be
invoked is used (see below):

    $ phylo_pipeline <STEP> [OPTIONS] <MAKEFILE>


## 3.2 Samples
The phylogenetic pipeline expects a number of samples
to be specified. Each sample has a name, a gender, and a genotyping
method:

    Samples:
      <GROUP>:
        SAMPLE_NAME:
          Gender: ...
          Genotyping Method: ...

Gender is required, and is used to filter SNPs at homozygous sex
chromsomes (e.g. chrX and chrY for male humans). Any names may be
used, and can simply be set to e.g. 'NA' in case this feature is not
used.

The genotyping method is either "SAMTools" for the default genotyping
procedure using samtools mpileupe | bcftools view, or "Random
Sampling" to sample one random nucleotide in the pileup at each
position. This key may be left out to use the default (SAMTools)
method.

Groups are optional, and may be used either for the sake of the
reader, or to specify a group of samples in lists of samples,
e.g. when excluding samples form a subsequent step, when filtering
singletons, or when rooting phylogenetic trees (see below)

For a given sample with name S, and a prefix with name P, the pipeline
will expect files to be located at ./data/samples/S.P.bam, or at
./data/samples/S.P.realigned.bam if the "Realigned" option is enabled
(see below).


## 3.3 Regions of interest
Analysis is carred out for a set of
"Regions of Interest", which is defined a set of named regions
specified using BED files:

    RegionsOfInterest:
      NAME:
        Prefix: NAME_OF_PREFIX
        Realigned: yes/no
	    ProteinCoding: yes/no
        IncludeIndels: yes/no

The options 'ProteinCoding' and 'IncludeIndels' takes values 'yes' and
'no' (without quotation marks), and determines the behavior when
calling indels. If 'IncludeIndels' is set to yes, indels are included
in the consensus sequence, and if 'ProteinCoding' is set to yes, only
indels that are a multiple of 3bp long are included.

The name and the prefix determines the location of the expected BED
file and the FASTA file for the prefix: For a region of interest named
R, and a prefix named P, the pipeline will expect the BED file to be
located at ./data/regions/P.R.bed. The prefix file is expected to be
located at ./data/prefixes/P.fasta


## 3.4 Genotyping

Genotyping is done either by random sampling of positions, or by
building a pileup using samtools and calling SNPs / indels using
bcftools. The command used for full genotyping is similar to the
following command:

    $ samtools mpileup [OPTIONS] | bcftools view [OPTIONS] -

In addition, SNPs / indels are filtered using the script 'vcf_filter',
which is included with the pipeline. This script implements the
filteres found in "vcfutils.pl varFilter", with some additions.

Options for either method, including for both "samtools mpileup" and
the "bcftools view" command is set using the **Genotyping** section of
the makefile, and may be set for all regions of interest (default
behavior) or for each set of regions of interest:

    Genotyping:
      Defaults: ...

The 'Defaults' key specifies that the options given here apply to all
regions of interest; in addition to this key, the name of each set of
regions of interest may be used, to set specific values for one set of
regions vs. another set. Thus, assuming regions of interest 'ROI\_a'
and 'ROI\_b', options may be set as follows:

    Genotyping:
      Defaults: ...
      ROI_a: ...
      ROI_b: ...

For each set of regions of interest named ROI, the final settings are
derived by first taking the Defaults, and then overwriting values
using the value taken from the ROI section (if one such exists). The
following shows how to change values in Defaults for a single ROI:

    Genotyping:
      Defaults:
        --switch: value_a
      ROI_N:
        --switch: value_b

In the above, all ROI except "ROI\_N" will use the switch with
'value\_a', while "ROI\_N" will use 'value\_b'. Executing the
'genotyping' step is described in 3.7.

Finally, note the "Padding" option; this option specifies a number of
bases to include around each interval in a set of regions of
interest. The purpose of this padding is to allow filtering of SNPs
based on the distance from indels, in the case where the indels are
outside the intervals themselves.


## 3.5 Multiple sequence alignment

Multiple sequence alignment (MSA) is currently carried out using
MAFFT, if enabled. Note that it is still nessesary to run the MSA
command (see below), even if the multiple sequence alignment itself is
disabled (for example in the case where indels are not called in the
genotyping step). This is because the MSA step is responsible for
generating both the unaligned multi-FASTA files, and the aligned
multi-FASTA files. It is nessesary to run the 'genotyping' step prior
to running the MSA step (see above).

It is possible to select among the various MAFFT algorithms using the
"Algorithm" key, and additionally to specify command-line options for
the selected algorithm:

    MultipleSequenceAlignment:
      Defaults:
        Enabled: yes

        MAFFT:
          Algorithm: G-INS-i
          --maxiterate: 1000

Currently supported algorithms are as follows (as described on the MAFFT
[website](http://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html)):

 * mafft - The basic program (mafft)
 * auto - Equivalent to command 'mafft --auto'
 * fft-ns-1 - Equivalent to the command 'fftns --retree 1'
 * fft-ns-2 - Equivalent to the command 'fftns'
 * fft-ns-i - Equivalent to the command 'fftnsi'
 * nw-ns-i - Equivalent to the command 'nwnsi'
 * l-ins-i - Equivalent to the command 'linsi'
 * e-ins-i - Equivalent to the command 'einsi'
 * g-ins-i - Equivalent to the command 'ginsi'

Command line options are specified as key / value pairs, as shown
above for the --maxiterate option, in the same manner that options are
specified for the genotyping section. Similarly, options may be
specified for all regions of interest ("Defaults"), or using the name
of a set of regions of interest, in order to set options for only that
set of regions.


## 3.6 Phylogenetic inference

Maximum likelyhood Phylogenetic inference is carried out using the
ExaML program. A phylogeny consists of a named (subsets of) one or
more sets of regions of interest, with individual regions partitioned
according to some scheme, and rooted on the midpoint of the tree or
one or more taxa:

    PhylogeneticInference:
      PHYLOGENY_NAME:
        ExcludeSamples:
		  ...

        RootTreesOn:
		  ...

        PerGeneTrees: yes/no

        RegionsOfInterest:
           REGIONS_NAME:
             Partitions: "111"
             SubsetRegions: SUBSET_NAME

        ExaML:
          Replicates: 1
          Bootstraps: 100
          Model: GAMMA


A phylogeny may exclude any number of samples specified in the Samples
region, by listing them under the ExcludeSamples. Furthermore, if
groups have been specified for samples (e.g. "<name>"), then these may
be used as a short-hand for multiple samples, by using the name of the
group including the angle-brackets ("<name>").

Rooting is determined using the RootTreesOn options; if this option is
not set, then the resulting trees are rooted on the midpoint of the
tree, otherwise it is rooted on the clade containing all the given
taxa. If the taxa does not form a monophyletic clade, then rooting is
done on the monophyletic clade containing the given taxa.

If PerGeneTrees is set to yes, a tree is generated for every named
feature in the regions of interest (e.g. genes), otherwise a
super-matrix is created based on all features in all the regions of
interest specified for the current phylogeny.

Each phylogeny may include one or more sets of regions of interest,
specified under the "RegionsOfInterest", using the same names as those
specified under the Project section. Each feature in a set of regions
of interest may be partitioned according to position specific
scheme. These are specified using a string of numbers (0-9), which is
then applied across the selected sequences to determine the model for
each position. For example, for the scheme "abc" and a given
nucleotide sequence, models are applied as follows:

    AAGTAACTTCACCGTTGTGA
	abcabcabcabcabcabcabcab

Thus, the default partitioning scheme ("111") will use the same model
for all positions, and is equivalent to the schemes "1", "11", "1111",
etc. Similarly, a per-codon-position scheme may be accomplished using
"123" or a similar string. In addition to numbers, the character 'X' may
be used to exclude specific positions in an alignment. E.g. to exclude
the third position in codons, use a string like "11X". Alternatively,
Partitions may be set to 'no' to disable per-feature partitions;
instead a single partition is used per set of regions of interest.

The options in the ExaML section specifies the number of bootstrap
trees to generate from the original supermatrix, the number of
phylogenetic inferences to carry out on the original supermatrix
(replicate), and the model used (c.f. the ExaML documentation).

The name (PHYLOGENY_NAME) is used to determine the location of the
resulting files, by default ./results/TITLE/phylogenies/NAME/. If
per-gene trees are generated, an addition two folders are used, namely
the name of the regions of interest, and the name of the gene /
feature.

For each phylogeny, the following files are generated:

* alignments.partitions - List of partitions used when running ExaML;
  the "reduced" file contains the same list of partitions, after empty
  columns (no called bases) have been excluded.
* alignments.phy - Super-matrix used in conjunction with the list of
  partitions when calling ExaML; the "reduced" file contains the same
  matrix, but with empty columns (no bases called) excluded.
* alignments.reduced.binary - The reduced supermatrix / partitions in
  the binary format used by ExaML.
* bootstraps.newick - List of bootstrap trees in Newick format, rooted
  as specified in the makefile.
* replicates.newick - List of phylogenies inferred from the full
  super-matrix, rooted as specified in the makefile.
* replicates.support.newick - List of phylogenies inferred from the
  full super-matrix, with support values calculated using the
  bootstrap trees, and rooted as specified in the makefile.


## 3.7 Executing the pipeline

The phylogenetic pipeline is excuted similarly to the BAM pipeline,
except that a command is provided for each step ('genotyping', 'msa',
and 'phylogeny'):

    $ phylo_pipeline <COMMAND> [OPTIONS] <MAKEFILE>

Thus, to execute the genotyping step, the following command is used:

    $ phylo_pipeline genotyping [OPTIONS] <MAKEFILE>

In addition, it is possible to run multiple steps by joining these
with the plus-symbol. To run both the 'genotyping' and 'msa' step at
the same time, use the following command:

    $ phylo_pipeline genotyping+msa  [OPTIONS] <MAKEFILE>


---------------------
# 4. Troubleshooting

TODO


---------------------
# 5. Contact

For questions, bug reports, and suggestions, please contact Mikkel Schubert at [MSchubert@snm.ku.dk](mailto:MSchubert@snm.ku.dk)
