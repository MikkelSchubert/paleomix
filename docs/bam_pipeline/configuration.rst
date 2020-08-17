.. highlight:: ini
.. _bam_configuration:


Configuring the BAM pipeline
============================

The BAM pipeline exposes a number options, including the maximum number of threads used, and the maximum number of threads used for individual programs, the location of JAR files, and more. These may be set using the corresponding command-line options (e.g. `--max-threads`). However, it is also possible to set default values for such options. This is accomplished by writing options in ~/.paleomix/bam_pipeline.ini::

    max-threads = 16
    log-level = warning
    jar-root = /home/username/install/jar_root
    bwa-max-threads = 1
    temp-root = /tmp/username/bam_pipeline
    jre-options =
    bowtie2-max-threads = 1

.. note::
    Options in the configuration file correspond directly to command-line options for the BAM pipeline, with leading dashes (--) are removed; for example, the command-line option `--max-threads` becomes `max-threads` in the configuration file, as shown above.

These values will be used by the pipeline, unless the corresponding option is also supplied on the command-line. I.e. if `max_threads` is set to 4 in the `bam_pipeline.ini` file, but the pipeline is run using `paleomix bam_pipeline --max-threads 10`, then the max threads value is set to 10.

.. note::
    If no value is given for `--max-threads` in ini-file or on the command-line, then the maximum number of threads is set to the number of CPUs available for the current host.
