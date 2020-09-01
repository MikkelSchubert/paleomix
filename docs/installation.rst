.. highlight:: Bash
.. _installation:


Installation
============

The following instructions will install PALEOMIX for the current user, but does not include specific programs required by the pipelines. For pipeline specific instructions, refer to the requirements sections for the :ref:`BAM <bam_requirements>`, the :ref:`Phylogentic <phylo_requirements>`, and the :ref:`Zonkey <zonkey_requirements>` pipeline.

The recommended way of installing PALEOMIX is by use of the `pip`_ package manager for Python 3. If `pip` is not installed, then please consult the documentation for your operating system. For Debian based operating systems, `pip` may be installed as follows::

    $ sudo apt-get install python3-pip

In addition, some libraries used by PALEOMIX may require additional development files, namely those for zlib, libbz2, liblzma, and for Python 3::

    $ sudo apt-get install libz-dev libbz2-dev liblzma-dev python3-dev

Once all requirements have been installed, PALEOMIX may be installed using `pip`::

    $ python3 -m pip install paleomix

To verify that the installation was carried out correctly, run the command `paleomix`::

    $ paleomix
    PALEOMIX - pipelines and tools for NGS data analyses
    Version: 1.2.14

    ...

If you have not previously used `pip`, then you may need to add the `pip` bin folder to your PATH and restart your terminal before running the `paleomix` command::

    $ echo 'export PATH=~/.local/bin:$PATH' >> ~/.bashrc


Self-contained installation
---------------------------

The recommended method for installing PALEOMIX is using a virtual environment. Doing so
allows different versions of PALEOMIX to be installed simultaneously and ensures that PALEOMIX and its dependencies are not affected by the addition or removal of other python modules.

This installation method requires the `venv` module. On Debian based systems, this module must be installed separately::

    $ sudo apt-get install python3-venv

Once `venv` is installed, creation of a virtual environment and installation of PALEOMIX may be carried out as shown here::

    $ python3 -m venv venv
    $ source ./venv/bin/activate
    $ (venv) pip install paleomix
    $ (venv) deactivate

Following successful completion of these commands, the PALEOMIX tools will be accessible in the `./venv/bin/` folder. However, as this folder also contains a copy of Python itself, it is not recommended to add it to your PATH. Instead, simply link the `paleomix` commands to a folder in your PATH. This can be accomplished as follows::

    $ mkdir -p ~/.local/bin/
    $ ln -s ${PWD}/venv/bin/paleomix ~/.local/bin/

If ~/.local/bin is not already in your PATH, then it can be added as follows:

    $ echo 'export PATH=~/.local/bin:$PATH' >> ~/.bashrc


Upgrading an existing installation
----------------------------------

Upgrade an existing installation of PALEOMIX, installed using the methods described above, may also be accomplished using pip. To upgrade a regular installation, simply run `pip` install with the --upgrade option::

    $ pip install --upgrade paleomix

To upgrade an installation a self-contained installation, activate the environment before calling `pip`::

    $ source ./venv/bin/activate
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

    $ conda create -c bioconda -n paleomix python=2.7 pip adapterremoval=2.3.1 samtools=1.9 picard=2.22.9 bowtie2=2.3.5.1 bwa=0.7.17 mapdamage2=2.0.9 r-base=3.5.1 r-rcpp=1.0.4.6 r-rcppgsl=0.3.7 r-gam=1.16.1 r-inline=0.3.15

Alternatively, you can use the `environment.yaml` file contained in the PALEOMIX github repository::

    $ curl https://raw.githubusercontent.com/MikkelSchubert/paleomix/master/paleomix_environment.yaml
    $ conda env create -f paleomix_environment.yaml

.. note::
    The above command(s) currently only contain the dependencies for the BAM pipeline

You can now activate the paleomix environment with::

    $ conda activate paleomix

PALEOMIX is not within the dependencies list above, so we can install this
_within_ the environment as explained above::

    $ (paleomix) pip install --user paleomix

PALEOMIX requires the Picard JAR file in a specific place, we can symlink the versions in your conda environment into the correct place::

    $ (paleomix) mkdir -p ~/install/jar_root/
    $ (paleomix) ln -s /<path>/<to>/miniconda2/envs/paleomix/share/picard-2.22.9-0/picard.jar ~/install/jar_root/

.. note::
    If you're unsure what your PALEOMIX conda environment path is, you can see this by running `conda env list`.

Once completed, you can test the environment works correctly using the pipeline test commands described in :ref:`examples`.

To deactivate the paleomix environment, simply run::

    $ conda deactivate

If you ever need to remove the entire environment, run the following command::

    $ rm /<path>/<to>/miniconda2/envs/paleomix/
