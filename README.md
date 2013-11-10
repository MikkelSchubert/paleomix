BAM / Phylogenetic pipelines
==========================================
TODO

---------------------
# 1. Installation instructions
The following instructions will use ~/install as the installation
directory for the pipelines and required tools, but any directory may
be used. The two pipelines (BAM and phylogenetic) have their own
set of requirements, as well as a common set of requirements.

Note that while it is possible to install and run the pipelines on
MacOSX (tested on Mountain Lion), this is not recommended.

## 1.1. Common instructions
The following instructions apply to both pipelines included in the
repository. The pipelines require [Python 2.7](http://python.org), in
addition to a few modules. [Git](http://git-scm.com/)  is recommended
to ease upgrading of the pipeline (see below).

### 1.1.1. Obtaining the pypeline source-code
To install a copy of the pipeline, create a git clone as follows:

    $ mkdir -p ~/install
    $ cd ~/install
    $ git clone https://github.com/MikkelSchubert/pypeline.git

or download a ZIP archive containing the repository:

    $ mkdir -p ~/install
    $ cd ~/install
    $ wget https://github.com/MikkelSchubert/pypeline/archive/master.zip
    $ mv pypeline-master pypeline

### 1.1.2 Installing the pipeline
The least invasive installation is accomplished by updating your
PYTHONPATH and PATH variables as follows (for Bash users):

    $ echo "export PYTHONPATH=\$PYTHONPATH:~/install/pypeline" >> ~/.bashrc
    $ echo "export PATH=\$PATH:~/install/pypeline/bin" >> ~/.bashrc

Alternatively, you can make a local install using the included 'setup.py' script:

    $ cd ~/install/pypeline
    $ python setup.py install --user

Finally, the pipeline can be installed for all users using the same script:

    $ cd ~/install/pypeline
    $ sudo python setup.py install

Regardless of installation method, you may need to start a new Bash
session for these changes to take effect.


### 1.1.3. Installing required modules
The pipeline requires [Pysam v0.7.4+](https://code.google.com/p/pysam/):

    $ curl -O https://pysam.googlecode.com/files/pysam-0.7.5.tar.gz
    $ tar xvzf pysam-0.7.5.tar.gz
    $ cd pysam-0.7.5
    $ python setup.py install --user


## 1.2 Installing programs required by the BAM pipeline

The BAM pipeline requires the following programs:

 * [SAMTools](http://samtools.sourceforge.net) v0.1.18+,
 * [AdapterRemoval](http://code.google.com/p/adapterremoval/) v1.5+
 (please cite Lindgreen 2013).
 * [BWA](http://bio-bwa.sourceforge.net/) v0.5.9+, or 0.6.2+, or v0.7.5+
 * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) 2.1.0+

Only one of BWA and Bowtie2 needs to be installed, corresponding to
the program that the user enables in the makefiles (see below).

  * [mapDamage v2.0.2](http://ginolhac.github.io/mapDamage/)

Note that the GNU Scientific Library (GSL) and the R packages listed
in the mapDamage instructions are required for rescaling support;
inline, gam, Rcpp, RcppGSL and ggplot2 (>=0.9.2). Use the following
commands to verify that GSL / R packages have been correctly installed:

    $ mapDamage --check-R-packages
    $ gsl-config

In addition, the JARs for the following programs are required:
 * [Picard tools](http://picard.sourceforge.net/) v1.82+
 * [GATK](http://www.broadinstitute.org/gatk/) (any version)


The Picard Tools JARs (.jar) files are expected to be located in
~/install/jar_root/ by default, but this behaviour may be changed using
either the --jar-root command-line option, or the global configuration
file (see below).

The GATK JAR is only required if the user wishes to carry out local
realignment near indels (recommended), and is expected to be placed
in the same folder as the Picard Tools JARs.

**WARNING:** If alignments are to be carried out against the human nuclear
genome, chromosomes MUST be ordered by their number for GATK to work!
See the GATK
[FAQ](http://www.broadinstitute.org/gatk/guide/topic?name=faqs#id1204)
for more information.

## 1.3 Installing programs required by the phylogenetic pipeline

The following programs are required for various parts of the pipeline.

**Genotyping:**

 * [Samtools v1.18+](http://samtools.sourceforge.net/)
 * [BEDTools v2.16.0+] (https://code.google.com/p/bedtools/)
 * [Tabix v0.2.5](http://samtools.sourceforge.net/)

Both the ‘tabix’ and the ‘bgzip’ executable from the Tabix package must be installed.

**Multiple Sequence Alignment:**

* [MAFFT v7+](http://mafft.cbrc.jp/alignment/software/)

Note that the pipeline requires that the algorithm-specific MAFFT
commands (e.g. ‘mafft-ginsi’, ‘mafft-fftnsi’). These are automatically
created by the ‘make install’ command.

**Phylogenetic Inference:**

* [RAxML v7.3.2+](https://github.com/stamatak/standard-RAxML)
* [ExaML v1.0.5+, parser v1.0.2+] (https://github.com/stamatak/ExaML)

The pipeline expects a single-threaded binary named ‘raxmlHPC’ for
RAxML. The pipeline expects the ExaML binary to be named ‘examl’, and
the parser binary to be named ‘examlParser’. Compiling and running
ExaML requires MPI (e.g. Open-MPI; http://www.open-mpi.org/). Both
programs offer a variety of makefiles suited for different
server-architectures and use-cases. If in doubt, use the
Makefile.SSE3.gcc makefiles:

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
pipeline. Running these projects is recommended, in order to verify
that the pipeline and required applications have been correctly
installed:

*BAM pipeline:*

    $ cd ~/install/pypeline/examples/bam_pipeline/
    $ bam_pipeline run 000_makefile.yaml

*Phylogenetic pipeline:*

    $ cd ~/install/pypeline/examples/bam_pipeline/
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
        Library_2:
          ...
      Sample_2:
        ...
    Target_2:
      ...

### 2.1.3. Additional options for prefixes
Each prefix may be assigned a label (see the template for a list of
possible values), which is useful if these two are aligned
separately. If both the 'mitochondrial' and the 'nuclear' labels are
specified, the summary file that is generated will contain information
aggregating the information for each of these prefixes under the name
'endogenous'. Each label may be specified any number of times. For
example:

    Prefixes:
      hg19:
        Path:  prefixes/hg19.fa
        Label: nuclear
      chrM:
        Path:  prefixes/rCRS.fa
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
        Path:  prefixes/hg19.fa
        AreasOfInterest:
          Exome: prefixes/hg19.exome.bed

### 2.1.4. Alternative input files
By default, all files specified for lanes are trimmed using
AdapterRemoval. If reads have already been trimmed, these can be
specified as follows, using the same read-types as in the
'ExcludeReads' options. For example:

    ED209:
      ED209:
        ACGTTA:
          Lane01:
            Single: reads/trimmed.fa.gz # Path to SE reads
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
The 'phylo\_pipeline mkfile' command can be used to create a makefile
template, as with the 'bam\_pipeline mkfile' command (see 3.1). This
makefile is used to specify the samples, regions of interest (to be
analysed), and options for the various programs.

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
The phylogenetic pipeline expects a number of samples to be
specified. Each sample has a name, a gender, and a genotyping
method:

    Samples:
      <GROUP>:
        SAMPLE_NAME:
          Gender: ...
          Genotyping Method: ...

Gender is required, and is used to filter SNPs at homozygous
sex chromsomes (e.g. chrX and chrY for male humans). Any names may be
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
Analysis is carred out for a set of "Regions of Interest", which is
defined a set of named regions specified using BED files:

    RegionsOfInterest:
	  NAME:
	    Prefix: ...
		Realigned: ...
		ProteinCoding: ...
		IncludeIndels: ...

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



## 3.5 Multiple sequence analysis



## 3.6 Phylogenetic inference



## 3.7 Executing the pipeline


---------------------
# 4. Troubleshooting


---------------------
# Citations

 * Langmead B, Salzberg SL. "Fast gapped-read alignment with Bowtie 2" Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.
 * Li H, Durbin R. "Fast and accurate short read alignment with Burrows-Wheeler transform" Bioinformatics. 2009 Jul 15;25(14):1754-60. doi:10.1093/bioinformatics/btp324.
 * Li H, *et al*. "The Sequence Alignment/Map format and SAMtools" Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352.
 * Lindgreen S. "AdapterRemoval: easy cleaning of next-generation sequencing reads" BMC Res Notes. 2012 Jul 2;5:337. doi: 10.1186/1756-0500-5-337.
 * McKenna A, *et al*. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data" Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110.
