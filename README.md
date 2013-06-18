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

The pypeline and pypeline/bin folders need to be added to your PYTHONPATH and PATH variables respectively. This may be accomplished as follows (for bash users):

    $ echo "export PYTHONPATH=\$PYTHONPATH:~/install/pypeline" >> ~/.bashrc
    $ echo "export PATH=\$PATH:~/install/pypeline/bin" >> ~/.bashrc

You may need to re-login for these changes to take effect.


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
The BAM pipeline requires [SAMTools](http://samtools.sourceforge.net) v0.1.18+, [AdapterRemoval](http://code.google.com/p/adapterremoval/) v1.4+ (1.5+ recommended for new projects!), and [BWA](http://bio-bwa.sourceforge.net/) v0.5.9+ / v0.7.5+ or [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 2.0.0+.  Only the aligner (BWA/Bowtie2) that is to be used (selected in the individual makefiles) needs to be installed.

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
