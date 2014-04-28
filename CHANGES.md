=========
Changelog
=========

Current
=============
  * Add 'paleomix' command, which provides interface for the various tools
    included in the PALEOMIX pipeline; this reduces the number of executables
    exposed by the pipeline, and allows for prerequisite checks to be done in
    one place.
  * Reworking version checking; add checks for JRE version (1.6+), for GATK
    (to check that the JRE can run it), and improved error messages for
    unidentified and / or outdated versions, and reporting of version numbers
    and requirements.
  * Dispose of hsperfdata_* folders created by certain JREs when using a
    custom temporary directory, when running Picard tools.
  * Cleanup of error-message displayed if Pysam version is outdated.
  * Ensure that file-handles are closed in the main process before subprocess
    execution, to ensure that these recieve SIGPIPE upon broken pipes.
  * Fix manifest, ensuring that all files are included in source distribution.
  * Fix regression in coverage / depths, which would fail if invoked for
    specific regions of interest.


Version 1.0.0
=============
  * Switching to more traditional version-number tracking.
