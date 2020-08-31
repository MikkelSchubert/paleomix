.. highlight:: bash
.. _zonkey_requirements:


Software Requirements
=====================

In addition to the requirements listed for the PALEOMIX pipeline itself in the :ref:`installation` section, the Zonkey pipeline requires that other pieces of software be installed:

* RScript from the `R`_ package, v3.3.3.
* SmartPCA from the `EIGENSOFT`_ package, v13050+ [Patterson2006]_, [Price2006]_
* `ADMIXTURE`_ v1.23 [Alexander2009]_
* `PLINK`_ v1.07 [Chang2015]_
* `RAxML`_ v8.2.9 [Stamatakis2006]_
* `SAMTools`_ v1.3.1 [Li2009b]_
* `TreeMix`_ v1.12 [Pickrell2012]_

The following R packages are required in order to carry out the plotting:

* `RColorBrewer`_
* `ape`_ [Paradis2004]_
* `ggplot2`_ [Wickham2009]_
* `ggrepel`_
* `reshape2`_ [Wickham2007]_

The R packages may be installed using the following commands::

    $ R
    > install.packages(c('RColorBrewer', 'ape', 'ggrepel', 'ggplot2', 'reshape2'))


Zonkey reference panel
----------------------

Running the Zonkey pipeline requires a reference panel containing the information needed for hybrid identification. A detailed description of the reference panel and instructions for where to download the latest version can be found in the :ref:`zonkey_panel` section.


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
