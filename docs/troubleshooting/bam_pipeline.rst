.. _troubleshooting_bam:

Troubleshooting the BAM Pipeline
================================

Troubleshooting BAM pipeline makefiles
--------------------------------------

**Path included multiple times in target**:

    This message is triggered if the same target includes one more more input files more than once::

        Error reading makefiles:
          MakefileError:
            Path included multiple times in target:
              - Record 1: Name: ExampleProject, Sample: Synthetic_Sample_1, Library: ACGATA, Barcode: Lane_1_001
              - Record 2: Name: ExampleProject, Sample: Synthetic_Sample_1, Library: ACGATA, Barcode: Lane_3_001
              - Canonical path 1: /home/username/temp/bam_example/data/ACGATA_L1_R1_01.fastq.gz
              - Canonical path 2: /home/username/temp/bam_example/data/ACGATA_L1_R2_01.fastq.gz

    This may be caused by using too broad wildcards, or simple mistakes. The message indicates the lane in which the files were included, as well as the "canonical" (i.e. following the resolution of symbolic links, etc.) path to each of the files. To resolve this issue, ensure that each input file is only included once for a given target.


**Target name used multiple times**:

    If running multiple makefiles in the same folder, it is important that the names given to targets in each makefile are unique, as the pipeline will otherwise mixfiles between different projects (see the section :ref:`bam_filestructure` for more information). The PALEOMIX pipeline attempts to detect this, and prevents the pipeline from running in this case::

        Error reading makefiles:
          MakefileError:
            Target name 'ExampleProject' used multiple times; output files would be clobbered!

**OutOfMemoryException (Picard Tools):**

    By default, the BAM pipeline will limit the amount of heap-space used by Java programs to 4GB (on 64-bit systems, JVM defaults are used on 32-bit systems), which may prove insufficient in some instances. This will result in the failing program terminating with a stacktrace, such as the following::

        Exception in thread "main" java.lang.OutOfMemoryError
        at net.sf.samtools.util.SortingLongCollection.<init>(SortingLongCollection.java:101)
        at net.sf.picard.sam.MarkDuplicates.generateDuplicateIndexes(MarkDuplicates.java:443)
        at net.sf.picard.sam.MarkDuplicates.doWork(MarkDuplicates.java:115)
        at net.sf.picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:158)
        at net.sf.picard.sam.MarkDuplicates.main(MarkDuplicates.java:97)


    To resolve this issue, increase the maximum amount of heap-space used using the "--jre-option" command-line option; this permits the passing of options to the Java Runtime Environment (JRE). For example, to increase the maximum to 8gb, run the BAM pipeline as follows::

        $ paleomix bam run --jre-option=-Xmx8g [...]


Troubleshooting AdapterRemoval
------------------------------

The AdapterRemoval task will attempt to determine verify the quality-offset specified in the makefile; if the contents of the file does not match the expected offset (i.e. contains quality scores that fall outside the range expected with that offset http://en.wikipedia.org/wiki/FASTQ_format\#Encoding, the task will be aborted.

**Incorrect quality offsets specified in makefile**:

    In case where the sequence data can be determined to contain FASTQ records with a different quality offset than that specified in the makefile, the task will be aborted with the message corresponding to the following::

        <AdapterRM (SE): 'data/TGCTCA_L1_R1_02.fastq.gz' -> 'ExampleProject/reads/Synthetic_Sample_1/TGCTCA/Lane_1_002/reads.*'>: Error occurred running command:
          Error(s) running Node:
            Temporary directory: '/path/to/temp/folder'

          FASTQ file contains quality scores with wrong quality score offset (33); expected reads with quality score offset 64. Ensure that the 'QualityOffset' specified in the makefile corresponds to the input.
          Filename = data/TGCTCA_L1_R1_02.fastq.gz

    Please verify the format of the input file, and update the makefile to use the correct QualityOffset before starting the pipeline.


**Input file contains mixed FASTQ quality scores**:

    In case where the sequence data can be determined to contain FASTQ records with a different quality scores corresponding to the both of the possible offsets (for example both "!" and "a"), the task will be aborted with the message corresponding to the following example::

        <AdapterRM (SE): 'data/TGCTCA_L1_R1_02.fastq.gz' -> 'ExampleProject/reads/Synthetic_Sample_1/TGCTCA/Lane_1_002/reads.*'>: Error occurred running command:
            Error(s) running Node:
              Temporary directory: '/path/to/temp/folder'

            FASTQ file contains quality scores with both quality offsets (33 and 64); file may be unexpected format or corrupt. Please ensure that this file contains valid FASTQ reads from a single source.
            Filename = 'data/TGCTCA_L1_R1_02.fastq.gz'

    This error would suggest that the input-file contains a mix of FASTQ records from multiple sources, e.g. resulting from the concatenation of multiple sets of data. If so, make use of the original data, and ensure that the quality score offset set for each is set correctly.


**Input file does not contain quality scores**:

    If the input files does not contain any quality scores (e.g. due to malformed FASTQ records), the task will terminate, as these are required by the AdapterRemoval program. Please ensure that the input files are valid FASTQ files before proceeding.

    Input files in FASTA format / not in FASTQ format:

    If the input file can be determined to be in FASTA format, or otherwise be determined to not be in FASTQ format, the task will terminate with the following message::

        <AdapterRM (SE): 'data/TGCTCA_L1_R1_02.fastq.gz' -> 'ExampleProject/reads/Synthetic_Sample_1/TGCTCA/Lane_1_002/reads.*'>: Error occurred running command:
          Error(s) running Node:
            Temporary directory: '/path/to/temp/folder'

          Input file appears to be in FASTA format (header starts with '>', expected '@'), but only FASTQ files are supported.
          Filename = 'data/TGCTCA_L1_R1_02.fastq.gz'

    Note that the pipeline only supports FASTQ files as input for the trimming stage, and that these have to be either uncompressed, gzipped, or bzipped. Other compression schemes are not supported at this point in time.


Troubleshooting BWA
-------------------

**BWA prefix generated using different version of BWA / corrupt index**:

    Between versions 0.5 and 0.6, BWA changed the binary format used to store the index sequenced produced using the command "bwa index". Version 0.7 is compatible with indexes generated using v0.6. The pipeline will attempt to detect the case where the current version of BWA does not correspond to the version used to generate the index, and will terminate if that is the case.

    As the two formats contain both contain files with the same names, the two formats cannot co-exist in the same location. Thus to resolve this issue, either create a new index in a new location, and update the makefile to use that location, or delete the old index files (path/to/prefix.fasta.*), and re-index it by using the command "bwa index path/to/prefix.fasta", or by simply re-starting the pipeline.

    However, because the filenames used by v0.6+ is a subset of the filenames used by v0.5.x, it is possible to accidentally end up with a prefix that appears to be v0.5.x to the pipeline, but in fact contains a mix of v0.5.x and v0.6+ files. This situation, as well as corruption of the index, may result in the following errors:

    1. [bwt_restore_sa] SA-BWT inconsistency: seq_len is not the same

    2. [bns_restore_core] fail to open file './rCRS.fasta.nt.ann'

    3. Segmentation faults when running 'bwa aln'; these are reported as "SIGSEGV" in the file pipe.errors

    If this occurs, removing the old prefix files and generating a new index is advised (see above).


Troubleshooting validation of BAM files
---------------------------------------

**Both mates are marked as second / first of pair**:

    This error message may occur during validation of the final BAM, if the input files specified for different libraries contained duplicates reads (*not* PCR duplicate). In that case, the final BAM will contain multiple copies of the same data, thereby risking a significant bias in downstream analyses.

    The following demonstrates this problem, using a contrieved example based on the examples/bam_example project included with the pipeline::

        $ paleomix bam run makefile.yaml
        [...]
        <Validate BAM: 'ExampleProject.rCRS.bam'>: Error occurred running command:
          Error(s) running Node:
            Temporary directory: '/path/to/temp/folder'

          Error(s) running Node:
            Return-codes: [1]
            Temporary directory: '/path/to/temp/folder'

            <Command = ['java', '-server', '-Xmx4g',
                        '-Djava.io.tmpdir=/tmp/bam_pipeline/9a5beba9-1b24-4494-836e-62a85eb74bf3',
                        '-Djava.awt.headless=true', '-XX:+UseSerialGC', '-jar',
                        '/home/research/tools/opt/jar_root/ValidateSamFile.jar',
                        'I=ExampleProject.rCRS.bam',
                        'IGNORE=MATE_NOT_FOUND', 'IGNORE=INVALID_QUALITY_FORMAT']
             Status  = Exited with return-code 1
             STDOUT  = '/path/to/temp/folder/rCRS.validated'
             STDERR* = '/path/to/temp/folder/pipe_java_20885232.stderr'
             CWD     = '/home/temp/bam_example'>

    Picard's ValidateSamfile prints the error messages to STDOUT, the location of which is indicated above::

        $ cat '/tmp/bam_pipeline/9a5beba9-1b24-4494-836e-62a85eb74bf3/rCRS.validated'
        ERROR: Record 684, Read name Seq_101_1324_104_rv_0\2, Both mates are marked as second of pair
        ERROR: Record 6810, Read name Seq_1171_13884_131_fw_0\2, Both mates are marked as second of pair

    To identify the source of the problems, the problematic reads may be extracted from the BAM file::

        $ samtools view ExampleProject.rCRS.bam|grep -w "^Seq_101_1324_104_rv_0"
        Seq_101_1324_104_rv_0\2 131 NC_012920_1 1325 60 100M = 1325 -1 [...]
        Seq_101_1324_104_rv_0\2 131 NC_012920_1 1325 60 100M = 1325 1 [...]
        Seq_101_1324_104_rv_0\1 16 NC_012920_1 1327 37 51M2D49M * 0 0 [...]
        Seq_101_1324_104_rv_0\1 89 NC_012920_1 1327 60 51M2D49M * 0 0 [...]


    Note that both mate pairs are duplicated, with slight variations in the flags. The source of the reads may be determined using the "RG" tags (not shown here), which for files produced by the pipeline corresponds to the library names. Once these are known, the corresponding FASTQ files may be examined to determine the source of the duplicate reads. This problem should normally be detected early in the pipeline, as checks for the inclusion of duplicate data has been implemented (see below).

**Read ... found in multiple files**:

    In order to detect the presence of data that has been included multiple times, e.g. due to incorrect merging of data, the pipeline looks for alignments with identical names, sequences and quality scores. If such reads are found, the follow error is reported::

        <Detect Input Duplication: 15 files>: Error occurred running command:
          Read 'Seq_junk_682_0' found in multiple files:
            - 'ExampleProject/rCRS/Synthetic_Sample_1/ACGATA/Lane_1_002/paired.minQ0.bam'
            - 'ExampleProject/rCRS/Synthetic_Sample_1/ACGATA/Lane_1_001/paired.minQ0.bam'

           This indicates that the same data files have been included multiple times in the project. Please review the input files used in this project, to ensure that each set of data is included only once.

    The message given indicates which files (and hence which samples/libraries and lanes were affected, as described in section :ref:`bam_filestructure`). If only a single file is given, this suggests that the reads were also found in that one file.

    This problem may result from the accidental concatenation of files provided to the pipeline, or from multiple copies of the same files being included in the wildcards specified in the makefile. As including the same sequencing reads multiple times are bound to bias downstream analyses (if it does not cause validation failure, see sub-section above), this must be fixed before the pipeline is re-started.
