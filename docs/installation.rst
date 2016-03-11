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


Following succesful completion of these commands, the paleomix tools will be accessible in the ~/install/virtualenvs/paleomix/bin/ folder. However, as this folder also contains a copy of Python itself, it is not recommended to add it to your PATH. Instead, simply link the paleomix commands to a folder in your PATH. This can, for example, be accomplished as follows::

    $ mkdir ~/bin/
    $ echo 'export PATH=~/bin:$PATH' >> ~/.bashrc
    $ ln -s ~/install/virtualenvs/paleomix/bin/paleomix ~/bin/

PALEOMIX also includes a number of optional shortcuts which may be used in place of running 'paleomix <command>' (for example, the command 'bam_pipeline' is equivalent to running 'paleomix bam_pipeline')::

    $ ln -s ~/install/virtualenvs/paleomix/bin/bam_pipeline ~/bin/
    $ ln -s ~/install/virtualenvs/paleomix/bin/conv_gtf_to_bed ~/bin/
    $ ln -s ~/install/virtualenvs/paleomix/bin/phylo_pipeline ~/bin/
    $ ln -s ~/install/virtualenvs/paleomix/bin/bam_rmdup_collapsed ~/bin/
    $ ln -s ~/install/virtualenvs/paleomix/bin/trim_pipeline ~/bin/


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


.. warning::
    When upgrading to v1.2.x or later from version 1.1.x or an before, it is nessesary to perform a manual installation the first time. This is accomplished by downloading and unpacking the desired version of PALEOMIX from the list of releases, and then invoking setup.py. For example::

        $ wget https://github.com/MikkelSchubert/paleomix/archive/v1.2.3.tar.gz
        $ tar xvzf v1.2.3.tar.gz
        $ paleomix-1.2.3/
        # Either for the current user:
        $ python setup.py install --user
        # Or, for all users:
        $ sudo python setup.py install

    Once this has been done once, pip may be used to perform future upgrades as described above.


.. _pip: https://pip.pypa.io/en/stable/
.. _Pysam: https://github.com/pysam-developers/pysam/
.. _Python: http://www.python.org/
.. _virtualenv: https://virtualenv.readthedocs.org/en/latest/