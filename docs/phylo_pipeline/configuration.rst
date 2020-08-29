.. highlight:: ini
.. _phylo_configuration:


Configuring the phylogenetic pipeline
=====================================

The Phylo pipeline supports a number of command-line options (see `paleomix phylo --help`). These options may be set directly on the command-line (e.g. using `--max-threads 16`), but it is also possible to set default values for such options.

This is accomplished by writing options in `~/.paleomix/phylo_pipeline.ini`::

    max-threads = 16
    examl-max-threads = 7
    log-level = warning
    temp-root = /tmp/username/phylo_pipeline

Options in the configuration file correspond directly to command-line options for the Phylo pipeline, with leading dashes removed. For example, the command-line option `--max-threads` becomes `max-threads` in the configuration file.

Options specified on the command-line take precedence over those in the `phylo_pipeline.ini` file. For example, if `max-threads` is set to 4 in the `phylo_pipeline.ini` file, but the pipeline is run using `paleomix phylo --max-threads 10`, then the max threads value is set to 10.
