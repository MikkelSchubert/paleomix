.. highlight:: ini
.. _bam_configuration:


Configuring the BAM pipeline
============================

The BAM pipeline supports a number of command-line options (see `paleomix bam run --help`). These options may be set directly on the command-line (e.g. using `--max-threads 16`), but it is also possible to set default values for such options.

This is accomplished by writing options in `~/.paleomix/bam_pipeline.ini`, such as::

    max-threads = 16
    bowtie2-max-threads = 1
    bwa-max-threads = 1
    jar-root = /home/username/install/jar_root
    jre-option = -Xmx4g
    log-level = warning
    temp-root = /tmp/username/bam_pipeline

Options in the configuration file correspond directly to command-line options for the BAM pipeline, with leading dashes removed. For example, the command-line option `--max-threads` becomes `max-threads` in the configuration file.

Options specified on the command-line take precedence over those in the `bam_pipeline.ini` file. For example, if `max-threads` is set to 4 in the `bam_pipeline.ini` file, but the pipeline is run using `paleomix bam run --max-threads 10`, then the max threads value is set to 10.
