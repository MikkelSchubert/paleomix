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
