.. highlight:: Bash
.. _installation:


Installation
============

The following instructions will install PALEOMIX for the current user, but does not include specific programs required by the pipelines. For pipeline specific instructions, refer to the requirements sections for the :ref:`BAM <bam_requirements>`, the :ref:`Phylogentic <phylo_requirements>`, and the :ref:`Zonkey <zonkey_requirements>` pipeline.

The recommended way of installing PALEOMIX is by use of the `pip`_ package manager for Python 3. If pip is not installed, then please consult the documentation for your operating system. For Debian based operating systems, pip may be installed as follows::

    $ sudo apt-get install python3-pip

In addition, some libraries used by PALEOMIX may require additional development files, namely those for `zlib`, `libbz2`, `liblzma`, and for Python 3::

    $ sudo apt-get install libz-dev libbz2-dev liblzma-dev python3-dev

Once all requirements have been installed, PALEOMIX may be installed using `pip`::

    $ python3 -m pip install paleomix==1.3.10

To verify that the installation was carried out correctly, run the command `paleomix`::

    $ paleomix
    PALEOMIX - pipelines and tools for NGS data analyses
    Version: v1.3.10

    ...

If you have not previously used pip, then you may need to add the pip `bin` folder to your `PATH` and restart your terminal before running the `paleomix` command::

    $ echo 'export PATH=~/.local/bin:$PATH' >> ~/.bashrc


Self-contained installation
---------------------------

The recommended method for installing PALEOMIX is using a virtual environment. Doing so
allows different versions of PALEOMIX to be installed simultaneously and ensures that PALEOMIX and its dependencies are not affected by the addition or removal of other python modules.

This installation method requires the `venv` module. On Debian based systems, this module must be installed separately::

    $ sudo apt-get install python3-venv

Once `venv` is installed, creation of a virtual environment and installation of PALEOMIX may be carried out as shown here::

    $ python3 -m venv venv
    $ ./venv/bin/pip install paleomix==v1.3.10

Following successful completion of these commands, the `paleomix` executable will be accessible in the `./venv/bin/` folder. However, as this folder also contains a copy of Python itself, it is not recommended to add it to your `PATH`. Instead, simply link the `paleomix` executable to a folder in your `PATH`. This can be accomplished as follows::

    $ mkdir -p ~/.local/bin/
    $ ln -s ${PWD}/venv/bin/paleomix ~/.local/bin/

If ~/.local/bin is not already in your PATH, then it can be added as follows::

    $ echo 'export PATH=~/.local/bin:$PATH' >> ~/.bashrc


Upgrading an existing installation
----------------------------------

Upgrade an existing installation of PALEOMIX, installed using the methods described above, may also be accomplished using pip. To upgrade a regular installation, simply run `pip install` with the `--upgrade` option::

    $ pip install --upgrade paleomix

To upgrade an installation a self-contained installation, simply call the `pip` executable in that environment::

    $ ./venv/bin/pip install --upgrade paleomix


Conda installation
------------------

`Conda`_ can be used to automatically setup a self-contained environment that includes the software required by PALEOMIX.

To install `conda` and also set it up so it can use the `bioconda`_ bioinformatics repository, follow the instructions on the bioconda website `here`_.

Next, run the following commands to download the conda environment template for this release of PALEOMIX and to create a new conda environment named `paleomix` using that template::

    $ curl -fL https://github.com/MikkelSchubert/paleomix/releases/download/v1.3.9/paleomix_environment.yaml > paleomix_environment.yaml
    $ conda env create -n paleomix -f paleomix_environment.yaml

You can now activate the paleomix environment with::

    $ conda activate paleomix

PALEOMIX requires that the Picard JAR file can be found in a specific location, so we can symlink the versions in your conda environment into the correct place::

    $ (paleomix) mkdir -p ~/install/jar_root/
    $ (paleomix) ln -s ~/*conda*/envs/paleomix/share/picard-*/picard.jar ~/install/jar_root/

.. note::
    If you installed conda in a different location, then you can obtain the location of the `paleomix` environment by running `conda env list`.

Once completed, you can test the environment works correctly using the pipeline test commands described in :ref:`examples`.

To deactivate the paleomix environment, simply run::

    $ conda deactivate

If you ever need to remove the entire environment, run the following command::

    $ conda env remove -n paleomix


.. _bioconda: https://bioconda.github.io
.. _conda: https://docs.conda.io/projects/conda/en/latest/index.html
.. _here: https://bioconda.github.io/user/install.html#install-conda
.. _pip: https://pip.pypa.io/en/stable/
.. _Pysam: https://github.com/pysam-developers/pysam/
.. _Python: http://www.python.org/
