.. highlight:: Bash
.. _bam_requirements:


Software requirements
=====================

In addition to the requirements listed in the :ref:`installation` section, the BAM pipeline requires several other pieces of software. The version numbers indicates the oldest supported version of each program:

* `AdapterRemoval`_ v2.2.0 [Schubert2016]_
* `SAMTools`_ v1.3.1 [Li2009b]_
* `Picard Tools`_ v1.137

The Picard Tools JAR-file (`picard.jar`) is expected to be located in `~/install/jar_root` by default, but this behavior may be changed using either the `--jar-root` command-line option, or via the global configuration file (see section :ref:`bam_configuration`)::

    $ mkdir -p ~/install/jar_root
    $ wget -O ~/install/jar_root/picard.jar https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar

Running Picard requires a Jave Runtime Environment (i.e. the `java` command). Please refer to your distro's documentation for how to install a JRE.

Furthermore, one or both of the following sequence aligners must be installed:

* `Bowtie2`_ v2.3.0 [Langmead2012]_
* `BWA`_ v0.7.15 [Li2009a]_

mapDamage is required by default, but can be disabled on a per-project basis:

* `mapDamage`_ 2.2.1 [Jonsson2013]_

If mapDamage is used to perform rescaling of post-mortem DNA damage, then the GNU Scientific Library (GSL) and the R packages listed in the mapDamage installation instructions are required; these include `inline`, `gam`, `Rcpp`, `RcppGSL` and `ggplot2`. Use the following commands to verify that these packages have been correctly installed::

    $ gsl-config
    Usage: gsl-config [OPTION]
    ...

    $ mapDamage --check-R-packages
    All R packages are present


On Debian-based systems, most of these dependencies can be installed using the following command::

    $ sudo apt-get install adapterremoval samtools bowtie2 bwa mapdamage

Testing the pipeline
--------------------

An example project is included with the BAM pipeline, and it is recommended to run this project in order to verify that the pipeline and required applications have been correctly installed. See the :ref:`examples` section for a description of how to run this example project.

.. Note::
    The example project does not carry out rescaling using mapDamage by default. If you wish to test that the requirements for mapDamage rescaling feature have been installed correctly, then change the following line

    .. code-block:: yaml

        mapDamage: plot

    to

    .. code-block:: yaml

        mapDamage: rescale

In case of errors, please consult the :ref:`troubleshooting` section.


.. _AdapterRemoval: https://github.com/MikkelSchubert/adapterremoval
.. _Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/
.. _BWA: http://bio-bwa.sourceforge.net/
.. _mapDamage: http://ginolhac.github.io/mapDamage/
.. _SAMTools: https://samtools.github.io
.. _Picard Tools: http://broadinstitute.github.io/picard/