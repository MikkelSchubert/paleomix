.. highlight:: Bash
.. _phylo_requirements:


Software requirements
=====================

Depending on the parts of the Phylogenetic pipeline used, different programs are required. The following lists which programs are required for each pipeline, as well as the minimum version required:


Genotyping
----------

* `BEDTools`_ v2.16.0+ [Quinlan2010]_
* [SAMTools](http://samtools.sourceforge.net) v0.1.18+ [Li2009b]_
* `Tabix`_ v0.2.5

Both the 'tabix' and the 'bgzip' executable from the Tabix package must be installed.


Multiple Sequence Alignment
---------------------------

* `MAFFT`_ v7+ [Katoh2013]_

Note that the pipeline requires that the algorithm-specific MAFFT commands (e.g. 'mafft-ginsi', 'mafft-fftnsi'). These are automatically created by the 'make install' command.


Phylogenetic Inference
----------------------

* `RAxML`_ v7.3.2+ [Stamatakis2006]_
* `ExaML`_ v1.0.5+

The pipeline expects a single-threaded binary named 'raxmlHPC' for RAxML. The pipeline expects the ExaML binary to be named 'examl', and the parser binary to be named 'parse-examl'. Compiling and running ExaML requires an MPI implementation (e.g. `OpenMPI`_), even if ExaML is run single-threaded. On Debian and Debian-based distributions, this may be accomplished installing 'mpi-default-dev' and 'mpi-default-bin'.

Both programs offer a variety of makefiles suited for different server-architectures and use-cases. If in doubt, use the Makefile.SSE3.gcc makefiles, which are compatible with most modern systems::

    $ make -f Makefile.SSE3.gcc


Testing the pipeline
--------------------

An example project is included with the phylogenetic pipeline, and it is recommended to run this project in order to verify that the pipeline and required applications have been correctly installed. Note that running this example project requires that the BAM pipeline has already been setup, as described in section :ref:`bam_requirements`):

    $ cd ~/install/pypeline/examples/phylo_pipeline/alignment/
    $ ./setup.sh
    $ bam_pipeline run 000_makefile.yaml
    $ cd ../phylogeny
    $ phylo_pipeline genotype+msa+phylogeny 000_makefile.yaml

The resulting files are described in section :ref:`phylo_filestructure`.


.. _BEDTools: https://code.google.com/p/bedtools/
.. _Tabix: http://samtools.sourceforge.net/
.. _MAFFT: http://mafft.cbrc.jp/alignment/software/
.. _RAxML: https://github.com/stamatak/standard-RAxML
.. _EXaML: https://github.com/stamatak/ExaML
.. _OpenMPI: http://www.open-mpi.org/