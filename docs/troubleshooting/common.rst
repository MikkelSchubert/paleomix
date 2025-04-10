.. highlight:: Bash
.. _troubleshooting_common:

Troubleshooting general problems
================================


If a command fails while the pipeline is running (e.g. mapping, genotyping, validation of BAMs, etc.), the pipeline will print a message to the command-line and write a message to a log-file. The location of the log-file may be specified using the --log-file command-line option, but if --log-file is not specified, a time-stamped log-file is generated in the temporary folder specified using the --temp-root command-line option, and the location of this log-file is printed by the pipeline during execution::

    $ 2014-01-07 09:46:19 Pipeline; 1 failed, 202 done of 203 tasks:
    Log-file located at '/path/to/temp/folder/bam_pipeline.20140107_094554_00.log'
    [...]


Most error-messages will involve a message in the following form::

    <Validate BAM: 'ExampleProject.rCRS.bam'>:
      Error occurred running command:
        Error(s) running Node:
          Return-codes: [1]
          Temporary directory: '/path/to/temp/folder'

          <Command = ['java', '-server', '-Xmx4g',
                      '-Djava.io.tmpdir=/path/to/temp/folder',
                      '-Djava.awt.headless=true', '-XX:+UseSerialGC', '-jar',
                      '/home/username/install/jar_root/ValidateSamFile.jar',
                      'I=ExampleProject.rCRS.bam',
                      'IGNORE=MATE_NOT_FOUND',
                      'IGNORE=INVALID_QUALITY_FORMAT']

           Status  = Exited with return-code 1
           STDOUT  = '/path/to/temp/folder/rCRS.validated'
           STDERR* = '/path/to/temp/folder/pipe_java_4454836272.stderr'
           CWD     = '/path/to/project'>

The task that failed was the validation of the BAM 'ExampleProject.rCRS.bam' using Picard ValidateSamFile, which terminated with return-code 1. For each command involved in a given task ('node'), the command-line (as the list passed to 'Popen'http://docs.python.org/2.7/library/subprocess.html), return code, and the current working directory (CWD) is shown. In addition, STDOUT and STDERR are always either piped to files, or to a different command. In the example given, STDOUT is piped to the file 'rCRS.validated', while STDERR is piped to the file 'pipe_java_4454836272.stderr'. The asterisks in 'STDERR*' indicates that this filename was generated by the pipeline itself, and that this file is only kept if the command failed.

To determine the cause of the failure (indicated by the non-zero return-code), examine the output of each command involved in the node. Normally, messages relating to failures may be found in the STDERR file, but in some cases (and in this case) the cause is found in the STDOUT file::

    $ cat /path/to/temp/folder/rCRS.validated
    ERROR: Record 87, Read name [...], Both mates are marked as second of pair
    ERROR: Record 110, Read name [...], Both mates are marked as first of pair
    [...]


This particular error indicates that the same reads have been included multiple times in the makefile (see section [sub:Troubleshooting-BAM]). Normally it is necessary to consult the documentation of the specified program in order to determine the cause of the failure.

In addition, the pipeline performs a number of which during startup, which may result in the following issues being detected:

**Required file does not exist, and is not created by a node**:

    Before start, the pipeline check for the presence of all required files. Should one or more files be missing, and the missing file is NOT created by the pipeline itself, an error similar to the following will be raised::

        $ paleomix bam run makefile.yaml
        [...]
        Errors detected during graph construction (max 20 shown):
          Required file does not exist, and is not created by a node:
            Filename: prefix/rCRS.fasta
            Dependent node(s): [...]

    This typically happens if the Makefile contains typos, or if the required files have been moved since the last time the makefile was executed. To proceed, it is necessary to determine the current location of the files in question, and/or update the makefile.


**Required executables are missing**:

    Before starting to execute a makefile, the pipeline will check that the requisite programs are installed, and verify that the installed versions meet the minimum requirements. Should an executable be missing, an error similar to the following will be issued, and the pipeline will not run::

        $ paleomix bam run makefile.yaml
        [...]
        Errors detected during graph construction (max 20 shown):
          Required executables are missing: bwa

    In that case, please verify that all required programs are installed (see sections TODO) and ensure that these are accessible via the current user's PATH (i.e. can be executed on the command-line using just the executable name).


**Version requirement not met**:

    In addition to checking for the presence of required executables (including java JARs), version of a program is checked. Should the version of the program not be compatible with the pipeline (e.g. because it is too old), the following error is raised::

        $ paleomix bam run makefile.yaml
        [...]
        Version requirement not met for 'Picard CreateSequenceDictionary.jar';
        please refer to the PALEOMIX documentation for more information.

            Executable:    /Users/mischu/bin/bwa
            Call:          bwa
            Version:       v0.5.7.x
            Required:      v0.5.19.x or v0.5.110.x or v0.6.2.x or at least v0.7.9.x

    If so, please refer to the documentation for the pipeline in question, and install/update the program to the version required by the pipeline. Note that the executable MUST be accessible by the PATH variable. If multiple versions of a program is installed, the version required by the pipeline must be first, which may be verified by using the "which" command::

        $ which -a bwa
        /home/username/bin/bwa
        /usr/local/bin/bwa

**Java Runtime Environment outdated / UnsupportedClassVersionError**:

    If the version of the Java Runtime Environment (JRE) is too old, the pipeline may fail to run with the follow message::

        The version of the Java Runtime Environment on this
        system is too old; please check the the requirement
        for the program and upgrade your version of Java.

        See the documentation for more information.

    Alternatively, Java programs may fail with a message similar to the following, as reported in the pipe_*.stderr file (abbreviated)::

        Exception in thread "main" java.lang.UnsupportedClassVersionError: org/broadinstitute/sting/gatk/CommandLineGATK :
          Unsupported major.minor version 51.0 at [...]

    To solve this problem, you will need to upgrade your copy of Java.
