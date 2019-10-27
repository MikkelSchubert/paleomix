# Changelog

## [1.2.13.8] - 2019-10-27
### Changed
  - Zonkey now identifies nuclear chromosomes by size instead of name; this
    is done to better handle FASTAs downloaded from different sources


## [1.2.13.7] - 2019-10-15
### Fixed
  - Fixed handling of digit only chromosome names in Zonkey
  - Remove dashes from Zonkey MT genomes when running 'mito' command


## [1.2.13.6] - 2019-10-13
### Fixed
  - Handle .*miss files created by some versions of plink in Zonkey


## [1.2.13.5] - 2019-09-29
### Fixed
  - Ignore ValidateSamFile warning REF_SEQ_TOO_LONG_FOR_BAI warning when
    processing genomes with contigs too large for BAI index files.


## [1.2.13.4] - 2019-03-25
### Fixed
  - Improved detection of Picard versions in cases where 'java'
    outputs additional text.


## [1.2.13.3] - 2018-11-01
### Fixed
  - Fixed validation/read counting of pre-trimmed reads not including
    the mate 1 files of paired-end reads. This resulted in the
    'seq_retained_reads' count being half the expected value.


## [1.2.13.2] - 2018-04-22
### Fixed
  - Additional fixes to divisions by zeros in summary calculations.
  - Fixed 'empty file' message if FASTA file ends with empty sequence.
  - Renamed pre-trimmed FASTQ validation/statistics file, to avoid failure if
    an older run was resumed.


## [1.2.13.1] - 2018-03-25
### Fixed 
  - Fixed divisions by zero if empty files are listed as pre-trimmed reads.


## [1.2.13] - 2018-03-25
### Added
  - Added 'retable' command for pretty-printing whitespace separated data in
    the previously used by the BAM pipeline.
  - Basic statistics are collected for pre-trimmed reads in the BAM pipeline.

### Changed
  - BAM pipeline tables are now saved as tab separated columns. The old
    pretty-printed format may be produced by running the 'retable' tool on
    the resulting files.
  - Memory usage for the 'coverage' and 'depths' commands were reduced when
    using very big BED files.

### Fixed
  - Fixed input / output files not being listed in 'pipe.errors' files.
  - Use the same max open files limit for picard (ulimit -n minus headroom)
    when determining if the default should be changed and as the final value.
  - Removed explicit test for JRE version, which was failing on some (valid)
    runtimes. Java programs are still checked prior to running pipelines.
  - Fixed changes to recent version of Pysam breaking the alignment step
    in the BAM pipeline.
  - Fixed various test failures resulting in different environments.
  - Fixed validation of pre-trimmed FASTQ files in BAM pipelines stopping
    early if empty files are encountered.

### Removed
  - Removed automatic migrating of configuration files created PALEOMIX
    PALEOMIX prior to v1.2.0.
  - Previously deprecated support for existing BAM files was removed from
    the BAM pipeline.


## [1.2.12] - 2017-08-13
### Fixed
  - Fixed input / output files not being listed in 'pipe.errors' files.
  - Use the same max open files limit for picard (ulimit -n minus headroom)
    when determining if the default should be changed and as the final value.

### Added
  - The 'vcf_to_fasta' command now supports VCFs containing haploid genotype
    calls, courtesy of Graham Gower.

### Changed
  - Require Pysam version 0.10.0 or later.


## [1.2.11] - 2017-06-09
### Fixed
  - Fixed unhandled exception if a FASTA file for a prefix is missing in a
    BAM pipeline makefile.
  - Fixed the 'RescaleQualities' option not being respected for non-global
    options in BAM pipeline makefiles.


## [1.2.10] - 2017-05-29
### Added
  - Preliminary support for CSI indexed BAM files, required for genomes with
    chromosomes > 2^29 - 1 bp in size. Support is still missing in HTSJDK, so
    GATK cannot currently be used with such genomes. CSI indexing is enabled
    automatically when required.

### Fixed
  - Reference sequences placed in the current directory no longer cause the
    BAM pipeline to complain about non-writable directories.
  - The maximum number of temporary files used by picard will no longer be
    increased above the default value used by the picard tools.

### Changed
  - The 'Status' of processes terminated by the pipeline will now be reported
    as 'Automatically terminated by PALEOMIX'. This is to help differentiate
    between processes that failed or were killed by an external source, and
    processes that were cleaned up by the pipeline itself.
  - Pretty-printing of commands shown when commands fail have been revised
    to make it more readable, including explicit descriptions when output
    is piped from one process to another and vice versa.
  - Commands are now shown in a format more suitable for running on the
    command-line, instead of as a Python list, when a node fails. Pipes are
    still specified separately.
  - Improved error messages for missing programs during version checks, and for
    exceptions raised when calling Popen during version checks.
  - Strip MC tags from reads with unmapped mates during cleanup; this is
    required since Picard (v2.9.0) ValidateSamFile considers such tags invalid.


## [1.2.9] - 2017-05-01
### Fixed
  - Improved handling of BAM tags to prevent unintended type changes.
  - Fixed 'rmdup_collapsed' underreporting the number of duplicate reads (in
    the 'XP' tag), when duplicates with different CIGAR strings were processed.

### Changed
  - PCR duplicates detected for collapsed reads using 'rmdup\_collapsed' are
    now identified based on alignments that include clipped bases. This
    matches the behavior of the Picard 'MarkDuplicates' command.
  - Depending on work-load, 'rmdup\_collapsed' may now run up to twice as fast.


## [1.2.8] - 2017-04-28
### Added
  - Added FILTER entry for 'F' filter used in vcf\_filter. This corresponds to
    heterozygous sites where the allele frequency was not determined.
  - Added 'dupcheck' command. This command roughly corresponds to the
    DetectInputDuplication step that is part of the BAM pipeline, and attempts
    to identify duplicate data (not PCR duplicates), by locating reads mapped
    to the same position, with the same name, sequence, and quality scores.
  - Added link to sample data used in publication to the Zonkey documentation.

### Changed
  - Only letters, numbers, and '-', '_', and '.' are allowed in sample-names
    used in Zonkey, in order to prevent invalid filenames and certain programs
    breaking on whitespace. Trailing whitespace is stripped.
  - Show more verbose output when building Zonkey pipelines.
  - Picard tools version 1.137 or later is now required by the BAM pipeline.
    This is nessesary as newer BAM files (header version 1.5) would fail to
    validate when using earlier versions of Picard tools.

### Fixed
  - Fixed validation nodes failing on output paths without a directory.
  - Fixed possible uncaught exceptions when terminating cat commands used by
    FASTQ validation nodes resulting in loss of error messages.
  - Fixed makefile validation failing with an unhandled TypeError if unhashable
    types were found in unexpected locations. For example, a dict found where a
    subset of strings were allowed. These now result in a proper MakeFileError.
  - Fixed user options in the 'BWA' section of the BAM Pipeline makefiles not
    being correctly applied when using the 'mem' or the 'bwasw' algorithms.
  - Fixed some unit tests failing when the environment caused getlogin to fail.


## [1.2.7] - 2017-01-03
### Added
  - PALEOMIX now includes the 'Zonkey' pipeline, a pipeline for
    detecting equine F1 hybrids from archeological remains. Usage
    is described in the documentation.

### Changed
  - The wrongly named per-sample option 'Gender' in the phylogenetic
    pipeline makefile has been replaced with a 'Sex' option. This does not
    break backwards compatibility, and makefiles using the old name will still
    work correctly.
  - The 'RescaleQualities' option has been merged with the 'mapDamage' Feature
    in the BAM pipeline makefile. The 'mapDamage' feature now takes the options
    'plot', 'model', and 'rescale', allowing more fine-grained control.

### Fixed
  - Fixed the phylogenetic pipeline complaining about missing sample genders
    (now sex) if no regions of interest had been specified. The pipeline will
    now complain about there being no regions of interest, instead.
  - The 'random sampling' genotyper would misinterpret mapping qualities 10
    (encoded as '+') and 12 (encoded as '-') as indels, resulting in the
    genotyping failing. These mapping qualities are now correctly ignored.


## [1.2.6] - 2016-10-12
### Changed
  - PALEOMIX now uses the 'setproctitle' for better compatibility; installing
    / upgraing PALEOMIX using pip (or equivalent tools) should automatically
    install this dependency.

### Fixed
  - mapDamage plots should not require indexed BAMs; this fixed missing file
    errors for some makefile configurations.
  - Version check for java did now works correctly for OpenJDK JVMs.
  - Pressing 'l' or 'L' to list the currently running tasks now correctly
    reports the total runtime of the pipeline, rather than 0s.
  - Fixed broken version-check in setup.py breaking on versions of python
    older than than 2.7, preventing meaningful message (patch by beeso018).
  - The total runtime is now correctly reported when pressing the 'l' key
    during execution of a pipeline.
  - The logger will automatically create the output directory if this does
    not already exist; previously logged messages could cause the pipeline
    to fail, even if these were not in themselves fatal.
  - Executables required executables for version checks are now included in
    the prior checks for missing executables, to avoid version-checks failing
    due to missing executables.

### Added
  - PALEOMIX will attempt to automatically limit the per-process maximum
    number of file-handles used by when invoking Picard tools, in order
    to prevent failures due to exceeding the system limits (ulimit -n).


## [1.2.5] - 2015-06-06
### Changed
  - Improved information capture when a node raises an unexpected exception,
    mainly for nodes implementing their own 'run' function (not CommandNodes).
  - Improved printing of the state of output files when using the command-line
    option --list-output-files. Outdated files are now always listed as
    outdated, where previously these could be listed as 'Missing' if the task
    in question was queued to be run next.
  - Don't attempt to validate prefixes when running 'trim_pipeline'; note that
    the structure of the Prefix section the makefile still has to be valid.
  - Reverted commit normalizing the strand of unmapped reads.
  - The commands 'paleomix coverage' and 'paleomix depths' now accept records
    lacking read-group information by default; these are record as <NA> in the
    sample and library columns. It is further possible to ignore all read-group
    information using the --ignore-readgroups command-line option.
  - The 'bam_pipeline mkfile' command now does limited validation of input
    'SampleSheet.csv', prints generated targets sorted alphabetically, and
    automatically generates unique names for identically named lanes. Finally,
    the target template is not included automatically generating a makefile.
  - The 'coverage' and 'depth' commands are now capable of processing files
    containing reads with and without read-groups, without requiring the use
    of the --ignore-readgroups command-line option. Furthermore, reads for
    which the read-group is missing in the BAM header are treated as if no
    readgroup was specified for that read.
  - The 'coverage' and 'depth' command now checks that input BAM files are
    sorted during startup and while processing a file.
  - Normalized information printed by different progress UIs (--progress-ui),
    and included the maximum number of threads allowed.
  - Restructured CHANGELOG based on http://keepachangelog.com/

### Fixed
  - Fixed mislabeling of BWA nodes; all were labeled as 'SE'.
  - Terminate read duplication checks when reaching the trailing, unmapped
    reads; this fixes uncontrolled memory growth when an alignment produces a
    large number of unmapped reads.
  - Fixed the pipeline demanding the existence of files from lanes that had
    been entirely excluded due to ExcludeReads settings.
  - Fixed some tasks needlessly depending on BAM files being indexed (e.g.
    depth histograms of a single BAM), resulting in missing file errors for
    certain makefile configurations.
  - Fixed per-prefix scan for duplicate input data not being run if no BAMs
    were set to be generated in the makefile, i.e. if both 'RawBAM' and
    'RealignedBAM' was set to 'off'.

### Deprecated
  - Removed the BAM file from the bam_pipeline example, and added deprecation
    warning; support for including preexisting BAMs will be removed in a
    future version of PALEOMIX.


## [1.2.4] - 2015-03-14
### Added
  - Included PATH in 'pipe.errors' file, to assist debugging of failed nodes.

### Fixed
  - Fix regression causing 'fixmate' not to be run on paired-end reads. This
    would occasionally cause paired-end mapping to fail during validation.


## [1.2.3] - 2015-03-11
### Added
  - Added the ability to the pipelines to output the list of input files
    required for a given makefile, excluding any file built by the pipeline
    itself. Use the --list-input-files command-line option to view these.

### Changed
  - Updated 'bam_pipeline' makefile template; prefixes and targets are
    described more explicitly, and values for the prefix are commented out by
    default. The 'Label' option is no included in the template, as it is
    considered deprecated.
  - Allow the 'trim_pipeline' to be run on a makefile without any prefixes;
    this eases use of this pipeline in the case where a mapping is not wanted.
  - Improved handling of unmapped reads in 'paleomix cleanup'; additional
    flags (in particular 0x2; proper alignment) are now cleared if the mate is
    unmapped, and unmapped reads are always represented on the positive strand
    (clearing 0x4 and / or 0x20).


## [1.2.2] - 2015-03-10
### Added
  - Documented work-arounds for problem caused when upgrading an old version of
    PALEOMIX (< 1.2.0) by using 'pip' to install a newer version, in which all
    command-line aliases invoke the same tool.
  - Added expanded description of PALEOMIX to README file.
  - The tool 'paleomix vcf_filter' can now clear any existing value in the
    FILTER column, and only record the result of running the filters
    implemented by this tool. This behavior may be enabled by running
    vcf_filter with the command-line option '--reset-filter yes'.

### Changed
  - Improved parsing of 'depths' histograms when running the phylogenetic
    pipeline genotyping step with 'MaxDepth: auto'; mismatches between the
    sample name in the table and in the makefile now only cause a warning,
    allowing for the common case where files depths were manually recalculated
    (and --target was not set), or where files were renamed.
  - The tool 'paleomix rmdup_collapsed' now assumes that ALL single-end reads
    (flag 0x1 not set) are collapsed. This ensures that pre-collapsed reads
    used in the pipeline are correctly filtered. Furthermore, reads without
    quality scores will be filtered, but only selected as the unique
    representative for a set of potential duplicates if no reads have quality
    scores. In that case, a random read is selected among the candidates.

### Fixed
  - Fixed failure during mapping when using SAMTools v1.x.


## [1.2.1] - 2015-03-08
### Changed
  - Remove dependency on BEDTools from the Phylogenetic pipeline.
  - Change paleomix.__version__ to follow PEP 0396.

### Fixed
  - Stop 'phylo_pipeline makefile' from always printing help text.
  - Fixed bug causing the phylo_pipeline to throw exception if no additional
    command-line arguments were given.
  - Allow simulation of reads for phylogenetic pipeline example to be executed
    when PALEOMIX is run from a virtual environment.


## [1.2.0] - 2015-02-24
This is a major revision of PALEOMIX, mainly focused on reworking the internals
of the PALEOMIX framework, as well as cleaning up several warts in the BAM
pipeline. As a result, the default makefile has changed in a number of ways,
but backwards compatibility is still retained with older makefiles, with one
exception. Where previously the 'FilterUnmappedReads' would only be in effect
when 'MinQuality' was set to 0, this option is now independent of the
'MinQuality' option.

In addition, it is now possible to install PALEOMIX via Pypi, as described in
the (partially) updated documentation now hosted on ReadTheDocs.

### Changed
  - Initial version of updated documentation hosted on ReadTheDocs, to replace
    documentation currently hosted on the repository wiki.
  - mapDamage files and models are now only kept in the
    {Target}.{Prefix}.mapDamage folder to simplify the file-structure;
    consequently, re-scaling can be re-done with different parameters by
    re-running the model step in these folders.
  - Rework BWA backtrack mapping to be carried out in two steps; this requires
    saving the .sai files (and hence more disk-space used by intermediate
    files, which can be removed afterwards), but allows better control over
    thread and memory usage.
  - Validate paths in BAM makefiles, to ensure that these can be parsed, and
    that these do not contain keys other than '{Pair}'.
  - The mapping-quality filter in the BAM pipeline / 'cleanup' command now only
    applies to mapped reads; consequently, setting a non-zero mapq value, and
    setting 'FilterUnmappedReads' to 'no' will not result in unmapped reads
    being filtered.
  - Improved the cleanup of BAM records following mapping, to better ensure
    that the resulting records follow the recommendations in the SAM spec. with
    regards to what fields / flags are set.
  - Configuration files are now expected to be located in ~/.paleomix or
    /etc/paleomix rather than ~/.pypeline and /etc/pypeline. To ensure
    backwards compatibility, ~/.pypeline will be migrated when a pipeline is
    first run, and replaced with a symbolic link to the new location.
    Furthermore, files in /etc/pypeline are still read, but settings in
    /etc/paleomix take precedence.
  - When parsing GTF files with 'gtf_to_bed', use either the attribute
    'gene_type' or 'gene_biotype', defaulting to the value  'unknown_genetype'
    if neither attribute can be found; also support reading of gz / bz2 files.
  - The "ExcludeReads" section of the BAM Pipeline makefile is now a dictionary
    rather a list of strings. Furthermore, 'Singleton' reads are now considered
    seperately from 'Single'-end reads, and may be excluded independently of
    those. This does not break backwards compatibility, but as a consequence
    'Single' includes both single-end and singleton reads when using old
    makefiles.
  - Added command-line option --nth-sample to the 'vcf_to_fasta' command,
    allowing FASTA construction from multi-sample VCFs; furthermore, if no BED
    file is specified, the entire genotype is constructed assuming that the VCF
    header is present.
  - Modify the FASTA indexing node so that SAMTools v0.1.x and v1.x can be used
    (added workaround for missing feature in v1.x).
  - The "Features" section of the BAM Pipeline makefile is now a dictionary
    rather than a list of strings, and spaces have been removed from feature
    names. This does not break backwards compatibility.
  - EXaML v3.0+ is now required; the name of the examl parser executable is
    required to be 'parse-examl' (previously expected to be 'examlParser'),
    following the name used by EXaML v3.0+.
  - Pysam v0.8.3+ is now required.
  - AdapterRemoval v2.1.5+ is now required; it is now possible to provide a
    list of adapter sequences using --adapter-list, and to specify the number
    of threads uses by AdapterRemoval via the --adapterremoval-max-threads
    command-line option.
  - Renamed module from 'pypeline' to 'paleomix' to aviod conflicts.
  - Improved handling FASTQ paths containing wildcards in the BAM pipeline,
    including additional checks to catch unequal numbers of files for paired-
    end reads.
  - Switch to setuptools in preperation for PyPI registration.
  - Avoid seperate indexing of intermediate BAMs when possible, reducing the
    total number of steps required for typical runs.
  - Restructure tests, removing (mostly unused) node tests.
  - Reworked sub-command handling to enable migration to setup-tools, and
    improved the safety of invoking these from the pipeline itself.
  - The output of "trim_pipeline mkfile" now includes the section for
    AdapterRemoval, which was previously mistakenly omitted.
  - Increased the speed of the checks for duplicate input data (i.e. the same
    FASTQ record(s) included multiple times in one or more files) by ~4x.

### Added
  - Paleomix v1.2.0 is now available via Pypi ('pip install paleomix').
  - Added command 'paleomix ena', which is designed to ease the preparation of
    FASTQ reads previously recorded in a BAM pipeline makefile for submission
    to the European Nucleotide Archive; this command is current unstable, and
    not available by default (see comments in 'main.py').
  - Exposed 'bam_pipeline remap' command, which eases re-mapping the hits
    identified against one prefix against other prefixes.
  - Added validation of BED files supplied to the BAM pipeline, and expand
    validation of BED files supplied to the Phylogenetic pipeline, to catch
    some cases that may cause unexpected behavior or failure during runtime.
  - Support SAMTools v1.x in the BAM pipeline; note, however, that the
    phylogenetic pipeline still requires SAMTools v0.1.19, due to major changes
    to BCFTools 1.x, which is not yet supported.
  - Modified 'bam_cleanup' to support SAMTools 1.x; SAMTools v0.1.19 or v1.x+
    is henceforth required by this tool.
  - The gender 'NA' may now be used for samples for which no filtering of sex
    chromosomes is to be carried out, and defaults to an empty set of
    chromsomes unless explicitly overridden.
  - Pipeline examples are now available following installation via the commands
    "bam_pipeline example" and "phylo_pipeline example", which copy the example
    files to a folder specified by the user.
  - Added ability to specify the maximum number of threads used by GATK;
    currently only applicable for training of indel realigner.

### Fixed
  - Ensured that only a single header is generated when using multiple threads
    during genotyping, in order to avoid issues with programs unable to handle
    multiple headers.
  - Information / error messages are now more consistently logged to stderr, to
    better ensure that results printed to stdout are not mixed with those.
  - Fixed bug which could cause the data duplication detection to fail when
    unmapped reads were included.
  - Fixed default values not being shown for 'vcf_filter --help'.
  - Fix 'vcf_filter' when using pysam v0.8.4; would raise exception due to
    changes to the VCF record class.

### Removed
  - Removed the 'paleomix zip' command, as this is no longer needed thanks to
    built-in gzip / bzip2 support in AdapterRemoval v2.
  - Removed commandline options --allow-missing-input-files,
    --list-orphan-files, --target, and --list-targets.


## [1.1.1] - 2015-10-10

### Changed
  - Detect the presence of carriage-returns ('\r') in FASTA files used as
    prefixes; these cause issues with some tools, and files should be converted
    using e.g. 'dos2unix' first.

### Fixed
  - Minor fix to help-text displayed as part of running information.

### Deprecated
  - AdapterRemoval v1.x is now considered deprecated, and support will be
    dropped shortly. Please upgrade to v2.1 or later, which can be found at
    https://github.com/MikkelSchubert/adapterremoval

### Removed
  - Dropped support for Picard tools versions prior to 1.124; this was
    nessesitated Picard tools merging into a single jar for all commands. This
    jar (picard.jar) is expected to be located in the --jar-root folder.


## [1.1.0] - 2015-09-08
### Added
  - Check that regions of interest specified in PhylogeneticInference section
    corresponds to those specified earlier in the makefile.
  - Added the ability to automatically read MaxReadDepth values from
    depth-histograms generated by the BAM pipeline to the genotyping step.
  - Add support for BWA algorithms "bwasw" and "mem", which are recommended for
    longer sequencing reads. The default remains the "backtrack" algorithm.
  - Include list of filters in 'vcf_filter' output and renamed these to be
    compatible with GATK (using ':' instead of '=').
  - Support for genotyping entire BAM (once, and only once), even if only a set
    of regions are to be called; this is useful in the context of larger
    projects, and when multiple overlapping regions are to be genotyped.
  - Added validation of FASTA files for the BAM pipeline, in order to catch
    serveral types of errors that may lead to failure during mapping.
  - Added options to BAM / Phylo pipelines for writing Dot-file of the full
    dependency tree of a pipeline.
  - Added the ability to change the number of threads, and more, while the
    pipeline is running. Currently, already running tasks are not terminated if
    the maximum number of threads is decreased. Press 'h' during runtime to
    list commands.
  - Support for AdapterRemoval v2.
  - Allow the -Xmx option for Java to be overridden by the user.

### Changed
  - Prohibit whitespace and parentheses in prefix paths; these cause problems
    with Bowtie2, due to the wrapper script used by this program.
  - Allow "*" as the name for prefixes, when selecting prefixes by wildcards.
  - Rework genotyping step to improve performance when genotyping sparse
    regions (e.g. genes), and to allow transparent parallelization.
  - Require BWA 0.5.9, 0.5.10, 0.6.2, or 0.7.9+ for BWA backtrack; other
    versions have never been tested, or are known to contain bugs that result
    in invalid BAM files.
  - The memory limit it no longer increased for 32-bit JREs by default, as the
    value used by the pipeline exceeded the maxmimum for this architecture.
  - Improved verification of singleton-filtering settings in makefiles.
  - Reworked the 'sample_pileup' command, to reduce the memory usage for larger
    regions (e.g. entire chromosomes) by an order of magnitude. Also fixed some
    inconsistency in the calculation of distance to indels, resulting in some
    changes in results.
  - Changed 'gtf_to_bed' to group by the gene biotype, instead of the source.

### Fixed
  - Fixed a bug preventing new tasks from being started immediately after a
    task had failed; new tasks would only be started once a task had finished,
    or no running tasks were left.
  - Fixed MaxDepth calculation being limited to depths in the range 0 .. 200.
  - Added workaround for bug in Pysam, which caused parsing of some GTF files
    to fail if these contained unquoted values (e.g. "exon_number 2;").
  - Fixed bug causing some tasks to not be re-run if the input file changed.
  - Fixed off-by-one error for coverages near the end of regions / contigs.
  - Ensure that the correct 'paleomix' wrapper script is called when invoking
    the various other tools, even if this is not located in the current PATH.
  - Parse newer SAMTools / BCFTools version strings, so that a meaningful
    version check failure can be reported, as these versions are not supported
    yet due to missing functionality.
  - Fix potential deadlock in the genotyping tool, which could occur if either
    of the invoked commands failed to start or crashed / were killed during
    execution.
  - Fixed error in which summary files could not be generated if two (or more)
    prefixes using the same label contained contigs with overlapping names but
    different sizes.
  - Fixed problems calculating coverage, depths, and others, when when using a
    user-provided BED without a name column.
  - Improved termination of child-processes, when the pipeline is interrupted.

### Deprecated
  - The 'mkfile' command has been renamed to 'makefile' for both pipelines; the
    old command is still supported, but considered deprecated.

### Removed
  - Dropped support for the "verbose" terminal output due to excessive
    verbosity (yes, really). The new default is "running" (previously called
    "quiet"), which shows a list of currently running nodes at every update.


## [1.0.1] - 2014-04-30
### Added
  - Add 'paleomix' command, which provides interface for the various tools
    included in the PALEOMIX pipeline; this reduces the number of executables
    exposed by the pipeline, and allows for prerequisite checks to be done in
    one place.
  - Added warning if HomozygousContigs contains contigs not included in any of
    the prefixes specified in the makefile.

### Changed
  - Reworking version checking; add checks for JRE version (1.6+), for GATK
    (to check that the JRE can run it), and improved error messages for
    unidentified and / or outdated versions, and reporting of version numbers
    and requirements.
  - Dispose of hsperfdata_* folders created by certain JREs when using a
    custom temporary directory, when running Picard tools.
  - Cleanup of error-message displayed if Pysam version is outdated.
  - Ensure that file-handles are closed in the main process before subprocess
    execution, to ensure that these recieve SIGPIPE upon broken pipes.
  - Improvements to handling of implicit empty lists in makefiles; it is now
    no longer required to explicitly specify an empty list. Thus, the following
    is equivalent assuming that the pipeline expects a list:
      ExplicitEmptyList: []
      ImplicitEmptyList:
  - Tweak makefile templates; the phylo makefile now specifies Male/Female
    genders with chrM and chrX; for the BAM pipeline the ROIs sub-tree and
    Label is commented out by default, as these are optional.
  - Reduced start-up time for bigger pipelines.

### Fixed
  - Fix manifest, ensuring that all files are included in source distribution.
  - Fix regression in coverage / depths, which would fail if invoked for
    specific regions of interest.
  - Fix bug preventing Padding from being set to zero when genotyping.


## [1.0.0] - 2014-04-16
### Changed
  - Switching to more traditional version-number tracking.


[Unreleased]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.8...HEAD
[1.2.13.8]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.7...v1.2.13.8
[1.2.13.7]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.6...v1.2.13.7
[1.2.13.6]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.5...v1.2.13.6
[1.2.13.5]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.4...v1.2.13.5
[1.2.13.4]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.3...v1.2.13.4
[1.2.13.3]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.2...v1.2.13.3
[1.2.13.2]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13.1...v1.2.13.2
[1.2.13.1]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.13...v1.2.13.1
[1.2.13]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.12...v1.2.13
[1.2.12]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.11...v1.2.12
[1.2.11]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.10...v1.2.11
[1.2.10]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.9...v1.2.10
[1.2.9]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.8...v1.2.9
[1.2.8]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.7...v1.2.8
[1.2.7]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.6...v1.2.7
[1.2.6]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.5...v1.2.6
[1.2.5]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.4...v1.2.5
[1.2.4]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.3...v1.2.4
[1.2.3]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.2...v1.2.3
[1.2.2]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.1...v1.2.2
[1.2.1]: https://github.com/MikkelSchubert/paleomix/compare/v1.2.0...v1.2.1
[1.2.0]: https://github.com/MikkelSchubert/paleomix/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/MikkelSchubert/paleomix/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/MikkelSchubert/paleomix/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/MikkelSchubert/paleomix/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/MikkelSchubert/paleomix/compare/v1.0.0-RC...v1.0.0
