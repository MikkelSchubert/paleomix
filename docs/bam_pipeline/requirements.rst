.. highlight:: Bash
.. _bam_requirements:


Software requirements
=====================

In addition to the requirements listed in the ref:`installation` section, the BAM pipeline requires that a several other pieces of software be installed. The plus-sign following version numbers are used to indicate that versions newer than that version are also supported:

* `AdapterRemoval`_ v2.1+ [Lindgreen2012]_
* `SAMTools`_ v0.1.18+ [Li2009b]_
* `Picard Tools`_ v1.137+

The Picard Tools JAR-file (picard.jar) is expected to be located in ~/install/jar_root/ by default, but this behavior may be changed using either the --jar-root command-line option, or via the global configuration file (see section :ref:`bam_configuration`).

Furthermore, one or both of the following sequence aligners must be installed:

  * `Bowtie2`_ v2.1.0+ [Langmead2012]_
  * `BWA`_ v0.5.9+, v0.6.2, or v0.7.9+ [Li2009a]_

In addition, the following packages are used by default, but can be omitted if disabled during runtime:

* `mapDamage`_ 2.0.2+ [Jonsson2013]_
* `Genome Analysis ToolKit`_ [McKenna2010]_

If mapDamage is used to perform rescaling of post-mortem DNA damage, then the GNU Scientific Library (GSL) and the R packages listed in the mapDamage installation instructions are required; these include 'inline', 'gam', 'Rcpp', 'RcppGSL' and 'ggplot2' (>=0.9.2). Use the following commands to verify that these packages have been correctly installed::

    $ gsl-config
    Usage: gsl-config [OPTION]
    ...

    $ mapDamage --check-R-packages
    All R packages are present

The GATK JAR is only required if the user wishes to carry out local realignment near indels (recommended), and is expected to be placed in the same folder as the Picard Tools JAR (see above).

The example projects included in the PALEOMIX source distribution may be used to test that PALEOMIX and the BAM pipeline has been correctly installed. See the :ref:`examples` section for more information.

In case of errors, please consult the :ref:`troubleshooting` section.


Testing the pipeline
--------------------

An example project is included with the BAM pipeline, and it is recommended to run this project in order to verify that the pipeline and required applications have been correctly installed. See the :ref:`examples` section for a description of how to run this example project.


.. _AdapterRemoval: https://github.com/MikkelSchubert/adapterremoval
.. _Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/
.. _BWA: http://bio-bwa.sourceforge.net/
.. _mapDamage: http://ginolhac.github.io/mapDamage/
.. _Genome Analysis ToolKit: http://www.broadinstitute.org/gatk/
.. _SAMTools: https://samtools.github.io
.. _Picard Tools: http://broadinstitute.github.io/picard/