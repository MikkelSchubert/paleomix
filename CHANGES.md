=========
Changelog
=========

Current
=============
  * Dispose of hsperfdata_* folders created by certain JREs when using a
    custom temporary directory, when running Picard tools.
  * Cleanup of error-message displayed if Pysam version is outdated.
  * Reworking version checking; add checks for JRE version (1.6+), for GATK
    (to check that the JRE can run it), and improved error messages for
    unidentified and / or outdated versions, and reporting of version numbers
    and requirements.
  * Ensure that file-handles are closed in the main process before subprocess
    execution, to ensure that these recieve SIGPIPE upon broken pipes.
  * Fix manifest, ensuring that all files are included in source distribution.


Version 1.0.0
=============
  * Switching to more traditional version-number tracking.
