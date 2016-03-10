.. highlight:: Bash
.. _troubleshooting_install:

Throubleshooting the installation
=================================

**Pysam / Cython installation fails with "Python.h: No such file or directory" or "pyconfig.h: No such file or directory"**:

    Installation of Pysam and Cython requires that Python development files are installed. On Debian based distributions, for example, this may be accomplished by running the following command::

        $ sudo apt-get install python-dev


**Pysam installation fails with "zlib.h: No such file or directory"**:

    Installation of Pysam requires that "libz" development files are installed. On Debian based distributions, for example, this may be accomplished by running the following command::

        $ sudo apt-get install libz-dev


**Command not found when attempting to run 'paleomix'**:

    By default, the PALEOMIX executables ('paleomix', etc.) are installed in ~/.local/bin. You must ensure that this path is included in your PATH::

        $ export PATH=$PATH:~/.local/bin

    To automatically apply this setting on sub-sequent logins (assuming that you are using Bash), run the following command::

        $ echo "export PATH=\$PATH:~/.local/bin" >> ~/.bash_profile


**PALEOMIX command-line aliases invokes wrong tools**:

    When upgrading an old PALEOMIX installation (prior to v1.2.x) using pip, the existence of old files may result in all command-line aliases ('bam\_pipeline', 'phylo\_pipeline', 'bam\_rmdup\_collapsed', etc.) invoking the same command (typically 'phylo_pipeline')::

        $ bam_pipeline makefile.yaml
        Phylogeny Pipeline v1.2.1

        [...]

    This can be solved by removing these aliases, and then re-installing PALEOMIX using 'pip', shown here for a system-wide install::

        $ sudo rm -v /usr/local/bin/bam_pipeline /usr/local/bin/conv_gtf_to_bed /usr/local/bin/phylo_pipeline /usr/local/bin/bam_rmdup_collapsed /usr/local/bin/trim_pipeline
        $ sudo python setup.py install

    Alternatively, this may be resolved by downloading and manually installing PALEOMIX::

        $ wget https://github.com/MikkelSchubert/paleomix/archive/v1.2.1.tar.gz
        $ tar xvzf v1.2.1.tar.gz
        $ paleomix-1.2.1/
        # Either for the current user:
        $ python setup.py install --user
        # Or, for all users:
        $ sudo python setup.py install