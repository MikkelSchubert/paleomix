.. highlight:: bash
.. _zonkey_requirements:


Software Requirements
=====================

The Zonkey pipeline requires PALEOMIX version 1.2.7 or later. In addition to the requirements listed for the PALEOMIX pipeline itself in the :ref:`installation` section, the Zonkey pipeline requires that other pieces of software be installed:

* RScript from the `R`_ package, v3.1+.
* SmartPCA from the `EIGENSOFT`_ package, v13050+ [Patterson2006]_, [Price2006]_
* `ADMIXTURE`_ v1.23+ [Alexander2009]_
* `PLINK`_ v1.7+ [Chang2015]_
* `RAxML`_ v7.3.2+ [Stamatakis2006]_
* `SAMTools`_ v0.1.19+ [Li2009b]_
* `TreeMix`_ v1.12+ [Pickrell2012]_

The following R packages are required in order to carry out the plotting:

* `RColorBrewer`_
* `ape`_ [Paradis2004]_
* `ggplot2`_ [Wickham2009]_
* `ggrepel`_
* `reshape2`_ [Wickham2007]_

The R packages may be installed using the following commands::

    $ R
    > install.packages(c('RColorBrewer', 'ape', 'ggrepel', 'ggplot2', 'reshape2'))


Installing under OSX
--------------------

Installing the Zonkey pipeline under OSX poses several difficulties, mainly due to SmartPCA. In the follow, it is assumed that the `Brew package manager`_ has been installed, as this greatly simplifies the installation of other, required pieces of software.

Firstly, install software and libraries required to compile SmartPCA::

    $ brew install gcc
    $ brew install homebrew/dupes/lapack
    $ brew install homebrew/science/openblas

In each case, note down the values indicated for LDFLAGS, CFLAGS, CPPFLAGS, etc.

Next, download and unpack the `EIGENSOFT`_ software. The following has been tested on EIGENSOFT version 6.1.1 ('EIG6.1.1.tar.gz').

To build SmartPCA it may further be nessesary to remove the use of the 'real-time' library::

    $ sed -e's# -lrt##' Makefile > Makefile.no_rt

Once you have done this, you can build SmartPCA using the locally copied libraries::

    $ env CC="/usr/local/opt/gcc/bin/gcc-6" LDFLAGS="-L/usr/local/opt/openblas/lib/" CFLAGS="-flax-vector-conversions -I/usr/local/opt/lapack/include/" make -f Makefile.no_rt

The above worked on my installation, but you may need to correct the variables using the values provided by Brew, which you noted down after running the 'install' command. You may also need to change the location of GGC set in the CC variable.


Testing the pipeline
--------------------

An example project is included with the BAM pipeline, and it is recommended to run this project in order to verify that the pipeline and required applications have been correctly installed. See the :ref:`examples_zonkey` section for a description of how to run this example project.


.. _ADMIXTURE: https://www.genetics.ucla.edu/software/admixture/
.. _EIGENSOFT: http://www.hsph.harvard.edu/alkes-price/software/
.. _PLINK: https://www.cog-genomics.org/plink2
.. _R: http://www.r-base.org/
.. _RAxML: https://github.com/stamatak/standard-RAxML
.. _RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
.. _SAMTools: https://samtools.github.io
.. _TreeMix: http://pritchardlab.stanford.edu/software.html
.. _ape: https://cran.r-project.org/web/packages/ape/index.html
.. _ggrepel: https://cran.r-project.org/web/packages/ggrepel/index.html
.. _ggplot2: https://cran.r-project.org/web/packages/ggplot2/index.html
.. _reshape2: https://cran.r-project.org/web/packages/reshape2/index.html
.. _Brew package manager: http://www.brew.sh
