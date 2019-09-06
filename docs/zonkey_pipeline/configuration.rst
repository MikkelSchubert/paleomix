.. highlight:: ini
.. _zonkey_configuration:

Configuring the Zonkey pipeline
===============================

Unlike the :ref:`bam_pipeline` and the :ref:`phylo_pipeline`, the :ref:`zonkey_pipeline` does not make use of makefiles. However, the pipeline does expose a number options, including the maximum number of threads used, various program parameters, and more. These may be set using the corresponding command-line options (e.g. --max-threads to set the maximum number of threads used during runtime). However, it is also possible to set default values for such options, including on a per-host bases. This is accomplished by executing the following command, in order to generate a configuration file at ~/.paleomix/zonkey.ini:

.. code-block:: bash

    $ paleomix zonkey --write-config


The resulting file contains a list of options which can be overwritten::

	[Defaults]
	max_threads = 1
	log_level = warning
	treemix_k = 0
	admixture_replicates = 1
	ui_colors = on
	downsample_to = 1000000

These values will be used by the pipeline, unless the corresponding option is also supplied on the command-line. I.e. if "max_threads" is set to 4 in the "zonkey.ini" file, but the pipeline is run using "paleomix zonkey run --max-threads 10", then the max threads value is set to 10.

.. note::
    Options in the configuration file correspond directly to command-line options for the BAM pipeline, with two significant differences: The leading dashes (--) are removed and any remaining dashes are changed to underscores (_); as an example, the command-line option --max-threads becomes max\_threads in the configuration file, as shown above.

It is furthermore possible to set specific options depending on the current host-name. Assuming that the pipeline was run on multiple servers sharing a single home directory, one might set the maximum number of threads on a per-server basis as follows::

    [Defaults]
    max_threads = 32
    [BigServer]
    max_threads = 64
    [SmallServer]
    max_threads = 16


The names used (here "BigServer" and "SmallServer") should correspond to the hostname, i.e. the value returned by the "hostname" command:

.. code-block:: bash

    $ hostname
    BigServer

Any value set in the section matching the name of the current host will take precedence over the 'Defaults' section, but can still be overridden by specifying the same option on the command-line, as described above.
