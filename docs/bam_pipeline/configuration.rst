.. highlight:: ini
.. _bam_configuration:


Configuring the BAM pipeline
============================

The BAM pipeline exposes a number options, including the maximum number of threads used, and the maximum number of threads used for individual programs, the location of JAR files, and more. These may be set using the corresponding command-line options (e.g. --max-threads). However, it is also possible to set default values for such options, including on a per-host bases. This is accomplished by excuting the following command, in order to generate a configuration file at ~/.paleomix/bam_pipeline.ini:

.. code-block:: bash

    $ paleomix bam_pipeline --write-config


The resulting file contains a list of options which can be overwritten::

    [Defaults]
    max_threads = 16
    log_level = warning
    jar_root = /home/username/install/jar_root
    bwa_max_threads = 1
    progress_ui = quiet
    temp_root = /tmp/username/bam_pipeline
    jre_options =
    bowtie2_max_threads = 1
    ui_colors = on


These values will be used by the pipeline, unless the corresponding option is also supplied on the command-line. I.e. if "max_threads" is set to 4 in the "bam_pipeline.ini" file, but the pipeline is run using "paleomix bam_pipeline --max-threads 10", then the max threads value is set to 10.

.. note::
    If no value is given for --max-threads in ini-file or on the command-line, then the maximum number of threads is set to the number of CPUs available for the current host.

It is furthermore possible to set specific options depending on the current host-name. Assuming that the pipeline was run on multiple servers sharing a single home folder, one might set the maximum number of threads on a per-server basis as follows::

    [Defaults]
    ...
    [BigServer]
    max_threads = 64
    [SmallServer]
    max_threads = 16


The names used (here "BigServer" and "SmallServer") should correspond to the hostname, i.e. the value returned by the "hostname" command:

.. code-block:: bash

    $ hostname
    BigServer

Any value set in the section matching the name of the current host will take precedence over the 'Defaults' section, but can still be overridden by specifying the same option on the command-line.