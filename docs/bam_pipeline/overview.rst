Overview of analytical steps
============================

During a typical analyses, the BAM pipeline will proceed through the following steps for each sample:


1. Initial steps

    Each prefix (reference sequences in FASTA format) is indexed using `samtools faidx` and using the short-read aligner configured for the current project.

2. Preprocessing of reads

    Adapter sequences, low quality bases and ambiguous bases are trimmed; overlapping paired-end reads are merged, and short reads are filtered using AdapterRemoval [Schubert2016]_.

3. Mapping of reads

    1. Processed reads resulting from the adapter-trimming / read-collapsing step above are mapped using the chosen short-read aligner (BWA or Bowtie2). The resulting BAMs are tagged using the information specified in the makefile (sample, library, lane, etc.).

    2. The records of the resulting BAM are updated using `samtools fixmate` to ensure that PE reads contain the correct information about the mate read.

    3. The BAM is sorted using `samtools sort`, indexed using `samtools index`, and validated using Picard `ValidateSamFile`.

    4. Finally, the records are updated using `samtools calmd` to ensure consistent reporting of the number of mismatches relative to the reference genome (BAM tag 'NM').

4. Filtering of duplicates, recalculation (rescaling) of quality scores, and validation

    1. If enabled, PCR duplicates are filtered using Picard `MarkDuplicates` for SE and PE reads and using `paleomix rmdup_collapsed` for collapsed reads (see the :ref:`other_tools` section). PCR filtering is carried out per library.

    2. If mapDamage based rescaling of quality scores is, quality scores of bases that are potentially the result of *post-mortem* DNA damage are recalculated using a damage model built using mapDamage2.0 [Jonsson2013]_.

    3. The resulting BAMs are indexed and validated using Picard `ValidateSamFile`. Mapped reads at each position of the alignments are compared using the query name, sequence, and qualities. If a match is found, it is assumed to represent a duplication of input data (see :ref:`troubleshooting_bam`).

5. Generation of final BAMs

    Each BAM in the previous step is merged into a final BAM file.

6. Statistics

    1. If the `Summary` feature is enable, a single summary table is generated for each target. This table summarizes the input data in terms of the raw number of reads, the number of reads following filtering / collapsing, the fraction of reads mapped to each prefix, the fraction of reads filtered as duplicates, and more.

    2. Coverage statistics and depth histograms are calculated for the intermediate and final BAM files using `paleomix coverage` and `paleomix depths`, if enabled. Statistics are calculated genome-wide and for any regions of interest specified by the user.

    3. If mapDamage is enabled, mapDamage plots are generated; if modeling or rescaling is enabled, a model of the post-mortem DNA damage is also generated.
