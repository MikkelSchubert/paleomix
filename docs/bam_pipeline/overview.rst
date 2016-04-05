Overview of analytical steps
============================

During a typical analyses, the BAM pipeline will proceed through the following steps. Note that the exact order in which each step is carried out during execution is not necessarily as shown below, since the exact steps depend on the user settings, and since the pipeline will automatically run steps as soon as possible:


1. Initial steps

    1. Each prefix (reference sequences in FASTA format) is indexed using either "bwa index" or "bowtie-build", depending on the configuration used.

    2. Each prefix is indexed using "samtools faidx".

    3. A sequence dictionary is built for each prefix using Picard BuildSequenceDictionary.jar

2. Preprocessing of reads

    2. Adapter sequences, low quality bases and ambiguous bases are trimmed; overlapping paired-end reads are merged, and short reads are filtered using AdapterRemoval [Lindgreen2012]_.

3. Mapping of reads

    1. Processed reads resulting from the adapter-trimming / read-collapsing step above are mapped using the chosen aligner (BWA or Bowtie2). The resulting alignments are tagged using the information specified in the makefile (sample, library, lane, etc.).

    2. The records of the resulting BAM are updated using "samtools fixmate" to ensure that PE reads contain the correct information about the mate read).

    3. The BAM is sorted using "samtools sort", indexed using "samtools index" (if required based on the current configuration), and validated using Picard ValidateSamFile.jar.

    4. Finally, the records are updated using "samtools calmd" to ensure consistent reporting of the number of mismatches relative to the reference genome (BAM tag 'NM').

4. Processing of preexisting BAM files

    1. Any preexisting BAM files are re-tagged using Picard 'AddOrReplaceReadGroups.jar' to match the tagging of other reads processed by the pipeline.

    2. The resulting BAM is sorted, updated using "samtools calmd", indexed using "samtools index" (if required), and validated using Picard ValidateSamFile.jar.

5. Filtering of duplicates, rescaling of quality scores, and validation

    1. If enabled, PCR duplicates are filtered using Picard MarkDuplicates.jar (for SE and PE reads) and "paleomix rmdup_collapsed" (for collapsed reads; see the :ref:`other_tools` section). PCR filtering is carried out per library.

    2. If "Rescaling" is enabled, quality scores of bases that are potentially the result of *post-mortem* DNA damage are recalculated using mapDamage2.0 [Jonsson2013]_.

    3. The resulting BAMs are indexed and validated using Picard ValidateSamFile.jar. Mapped reads at each position of the alignments are compared using the query name, sequence, and qualities. If a match is found, it is assumed to represent a duplication of input data (see :ref:`troubleshooting_bam`).

6. Generation of final BAMs

    1. If the "Raw BAM" feature is enabled, each BAM in the previous step is merged into a final BAM file.

    2. If the "Realigned BAM" feature is enabled, each BAM generated in the previous step is merged, and GATK IndelRealigner is used to perform local realignment around indels, to improve downstream analyses. The resulting BAM is updated using "samtools calmd" as above.

7. Statistics

    1. If the "Summary" feature is enable, a single summary table is generated for each target. This table summarizes the input data in terms of the raw number of reads, the number of reads following filtering / collapsing, the fraction of reads mapped to each prefix, the fraction of reads filtered as duplicates, and more.

    2. Coverages statistics are calculated for the intermediate and final BAM files using "paleomix coverage", depending on makefile settings. Statistics are calculated genome-wide and for any regions of interest specified by the user.

    3. Depth histograms are calculated using "paleomix depths", similar to coverage statistics, these statistics are genome-wide and for any regions of interest specified by the user.

    4. If the "mapDamage" feature or "Rescaling" is enabled, mapDamage plots are generated; if rescaling is enabled, a model of the post-mortem DNA damage is also generated.

    5. If the "DuplicateHist" feature is enabled, histograms of PCR duplicates are estimated for each library, for use with the 'preseq' tool[Daley2013]_, to estimate the complexity of the libraries.
