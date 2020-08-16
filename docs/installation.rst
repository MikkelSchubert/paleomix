.. highlight:: Bash
.. _installation:


Installation
============

The following instructions will install PALEOMIX for the current user, but does not include specific programs required by the pipelines. For pipeline specific instructions, refer to the requirements sections for the :ref:`BAM <bam_requirements>`, the :ref:`Phylogentic <phylo_requirements>`, and the :ref:`Zonkey <zonkey_requirements>` pipeline. The recommended way of installing PALEOMIX is by use of the `pip`_ package manager for Python. If Pip is not installed, then please consult the documentation for your operating system.

In addition to the `pip`_ package manager for Python, the pipelines require `Python`_ 2.7, and `Pysam`_ v0.8.3+, which in turn requires both Python and libz development files (see the :ref:`troubleshooting_install` section). When installing PALEOMIX using pip, Pysam is automatically installed as well. However, note that installing Pysam requires the zlib and Python 2.7 development files. On Debian based distributions, these may be installed as follows:

    # apt-get install libz-dev python2.7-dev

.. warning::
  PALEOMIX has been developed for 64 bit systems, and has not been extensively tested on 32 bit systems!


Regular installation
--------------------

The following command will install PALEOMIX, and the Python modules required to run it, for the current user only::

    $ pip install --user paleomix

To perform a system-wide installation, simply remove the --user option, and run as root::

    $ sudo pip install paleomix

To verify that the installation was carried out correctly, run the command 'paleomix'::

    $ paleomix
    PALEOMIX - pipelines and tools for NGS data analyses.
    Version: v1.0.1

    Usage: paleomix <command> [options]
    [...]

If the command fails, then please refer to the :ref:`troubleshooting` section.


Self-contained installation
---------------------------

In some cases, it may be useful to make a self-contained installation of PALEOMIX, *e.g.* on shared servers. This is because Python modules that have been installed system-wide take precendence over user-installed modules (this is a limitation of Python itself), which may cause problems both with PALEOMIX itself, and with its Python dependencies.

This is accomplished using `virtualenv`_ for Python, which may be installed using `pip`_ as follows::

    $ pip install --user virtualenv

or (for a system-wide installation)::

    $ sudo pip install virtualenv


The follow example installs paleomix in a virtual environmental located in *~/install/virtualenvs/paleomix*, but any location may be used::

    $ virtualenv ~/install/virtualenvs/paleomix
    $ source ~/install/virtualenvs/paleomix/bin/activate
    $ (paleomix) pip install paleomix
    $ (paleomix) deactivate


Following successful completion of these commands, the paleomix tools will be accessible in the ~/install/virtualenvs/paleomix/bin/ folder. However, as this folder also contains a copy of Python itself, it is not recommended to add it to your PATH. Instead, simply link the paleomix commands to a folder in your PATH. This can, for example, be accomplished as follows::

    $ mkdir ~/bin/
    $ echo 'export PATH=~/bin:$PATH' >> ~/.bashrc
    $ ln -s ~/install/virtualenvs/paleomix/bin/paleomix ~/bin/


Upgrading an existing installation
----------------------------------

Upgrade an existing installation of PALEOMIX, installed using the methods described above, may also be accomplished using pip. To upgrade a regular installation, simply run pip install with the --upgrade option, for a user installation::

    $ pip install --user --upgrade paleomix

Or for a system-wide installation::

    $ sudo pip install --upgrade paleomix

To upgrade an installation a self-contained installation, simply activate the environment before proceeding::

    $ source ~/install/virtualenvs/paleomix/bin/activate
    $ (paleomix) pip install --upgrade paleomix
    $ (paleomix) deactivate

.. _pip: https://pip.pypa.io/en/stable/
.. _Pysam: https://github.com/pysam-developers/pysam/
.. _Python: http://www.python.org/
.. _virtualenv: https://virtualenv.readthedocs.org/en/latest/

Conda installation
-------------------

To have a completely contained environment that includes all software dependencies, you can create a [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment.

To install conda and also set it up so it can use the [bioconda](https://bioconda.github.io) bioinformatics tool repository, you can follow the instructions on the bioconda website [here](https://bioconda.github.io/user/install.html#install-conda).

Once set-up, you can create a conda environment using the following commands::

    $ conda create -c bioconda -n paleomix python=2.7 pip adapterremoval=2.3.1 samtools=1.9 picard=2.22.9 bowtie2=2.3.5.1 bwa=0.7.17 mapdamage2=2.0.9 gatk=3.8 r-base=3.5.1 r-rcpp=1.0.4.6 r-rcppgsl=0.3.7 r-gam=1.16.1 r-inline=0.3.15

Alternatively, you can use the `environment.yaml` file contained in the PALEOMIX github repository.

    $ curl https://raw.githubusercontent.com/MikkelSchubert/paleomix/master/paleomix_environment.yaml
    $ conda env create -f paleomix_environment.yaml

> Note the above command(s) currently only contain the dependencies for the bam_pipeline

You can now activate the paleomix environment with::

    $ conda activate paleomix

Paleomix is not within the dependencies list above, so we can install this
_within_ the environment as explained above::

    $ (paleomix) pip install --user paleomix

The bam_pipeline also needs older versions of GATK, which are now not maintained by the Broad Institute. We can download the JAR file from the Broad archive, and activate
it within the conda environment like so::

    $ (paleomix) wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
    $ (paleomix) gatk3-register GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

Paleomix requires the GATK and Picard JAR files in a specific place, we can symlink the versions in your conda environment into the correct place::

    $ (paleomix) mkdir -p /home/<YOUR_USER>/install/jar_root/
    $ (paleomix) ln -s /<path>/<to>/miniconda2/envs/paleomix/opt/gatk-3.8/GenomeAnalysisTK.jar /home/<user>/install/jar_root/
    $ (paleomix) ln -s /<path>/<to>/miniconda2/envs/paleomix/share/picard-2.22.9-0/picard.jar /home/<user>/install/jar_root/

> If you're unsure what your paleomix conda environment path is, you can see this by running `conda env list`.

Once completed, you can test the environment works correctly using the pipeline test commands described in :ref:`examples`.

To deactivate the paleomix environment, simply run::

    $ conda deactivate

If you ever need to remove the entire environment, run the following command::

    $ rm /<path>/<to>/miniconda2/envs/paleomix/
