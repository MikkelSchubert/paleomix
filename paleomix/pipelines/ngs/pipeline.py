import logging
import os

from paleomix.common.command import InputFile
from paleomix.common.fileutils import swap_ext
from paleomix.common.layout import Layout
from paleomix.nodes.adapterremoval import IdentifyAdaptersNode
from paleomix.nodes.bwa import BWAAlgorithmNode, BWAIndexNode, BWAMem2IndexNode
from paleomix.nodes.commands import FilterCollapsedBAMNode, FinalizeBAMNode
from paleomix.nodes.fastp import FastpNode
from paleomix.nodes.fastqc import FastQCNode
from paleomix.nodes.gatk import (
    ApplyBQSRNode,
    ApplyVQSRNode,
    BaseRecalibratorNode,
    CombineGVCFsNode,
    CreateSequenceDictionaryNode,
    FastqToSamNode,
    GatherVcfsNode,
    GenotypeGVCFs,
    HaplotypeCallerNode,
    SplitIntervalsNode,
    VariantRecalibratorNode,
)
from paleomix.nodes.multiqc import MultiQCNode
from paleomix.nodes.samtools import (
    BAMIndexNode,
    BAMStatsNode,
    FastaIndexNode,
    MarkDupNode,
    TabixIndexNode,
)
from paleomix.nodes.validation import ValidateFASTAFilesNode
from paleomix.pipelines.ngs.nodes import TranchesPlotsNode

########################################################################################
# Project layout

# File structure with leaf values representing labels used to refer a given path
_LAYOUT = {
    "alignments": {
        "{sample}.{genome}.bam": "aln_recal_bam",
        "{sample}.{genome}.junk.bam": "aln_split_failed_bam",
    },
    "cache": {
        "{sample}": {
            "alignments": {
                "{sample}.{genome}.{library}.{run}.{kind}.bam": "aln_run_bam",
                "{sample}.{genome}.{library}.{run}.{kind}.ValidateSamfile.log": "aln_run_validation_log",
                "{sample}.{genome}.{library}.rmdup.merged.bam": "aln_rmdup_merged_bam",
                "{sample}.{genome}.{library}.rmdup.paired.bam": "aln_rmdup_paired_bam",
                "{sample}.{genome}.{library}.rmdup.paired.metrics.txt": "aln_rmdup_paired_metrics",
                "{sample}.{genome}.good.bam": "aln_split_passed_bam",
                "{sample}.{genome}.good.ValidateSamfile.log": "aln_split_good_validation_log",
                "{sample}.{genome}.good.recalibration.table.txt": "aln_recal_training_table",
                "{sample}.{genome}.good.recalibration.table.log": "aln_recal_training_log",
                "{sample}.{genome}.good.recalibration.log": "aln_recal_log",
            },
            "reads": {
                "{sample}.{library}.{run}.paired_1.fastq.gz": "fastp_paired_1",
                "{sample}.{library}.{run}.paired_2.fastq.gz": "fastp_paired_2",
                "{sample}.{library}.{run}.fastq.gz": "fastp_merged",
                "{sample}.{library}.{run}.unpaired.fastq.gz": "fastp_unpaired",
                "{sample}.{library}.{run}.unpaired.bam": "fastp_unpaired_bam",
                "{sample}.{library}.{run}.failed.fastq.gz": "fastp_failed",
                "{sample}.{library}.{run}.failed.bam": "fastp_failed_bam",
            },
        },
        "haplotypes": {
            "{genome}.uncalibrated": {
                "{part}.g.vcf.gz": "gvcf_merged_part",
                "{part}.vcf.gz": "vcf_merged_part",
            },
            "{genome}.uncalibrated.g.vcf.gz": "gvcf_merged",
            "{genome}.uncalibrated.vcf.gz": "vcf_merged",
            "{genome}.recalibration.training.snp.vcf.gz": "vcf_recal_training_snp_vcf",
            "{genome}.recalibration.training.snp.vcf.gz.log": "vcf_recal_training_snp_vcf_log",
            "{genome}.recalibration.training.indel.vcf.gz": "vcf_recal_training_indel_vcf",
            "{genome}.recalibration.training.indel.vcf.gz.log": "vcf_recal_training_indel_vcf_log",
            "{genome}.recalibration.snp.vcf.gz": "vcf_recal_snp",
            "{genome}.recalibration.snp.log": "vcf_recal_snp_log",
            "{genome}.recalibration.snp.indel.log": "vcf_recal_snp_indel_log",
        },
    },
    "haplotypes": {
        "{sample}.{genome}.g.vcf.gz": "gvcf_per_sample",
        "{sample}.{genome}.g.vcf.gz.log": "gvcf_per_sample_log",
    },
    "genotypes": {
        "{genome}.vcf.gz": "vcf_recal_snp_indel",
    },
    "statistics": {
        "alignments_multiQC": "bam_multiqc_prefix",
        "alignments": {
            "{sample}.{genome}.mapping.json": "aln_split_statistics",
            "{sample}.{genome}.{method}.txt": "bam_stats",
            "fastqc": "bam_fastqc_dir",
        },
        "genotyping": {
            "{genome}.recalibration.snp.r": "vcf_recal_training_snp_r",
            "{genome}.recalibration.snp.tranches": "vcf_recal_training_snp_trances",
            "{genome}.recalibration.snp.tranches.extras": "vcf_recal_training_snp_trances_extras",
            "{genome}.recalibration.indel.r": "vcf_recal_training_indel_r",
            "{genome}.recalibration.indel.tranches": "vcf_recal_training_indel_trances",
        },
        "reads_pre_trimmed_multiQC": "stats_fastqc_multiqc_pre",
        "reads_trimming_multiQC": "stats_fastp_multiqc",
        "reads_post_trimmed_multiQC": "stats_fastqc_multiqc_post",
        "reads": {
            "fastp": {
                "{sample}": {
                    "{library}": {
                        "{run}.html": "fastp_report_html",
                        "{run}_fastp.json": "fastp_report_json",
                    },
                },
            },
            "pre_trimmed": {
                "{sample}": {
                    "{library}": {
                        "{run}": "stats_fastqc_pre",
                        "{run}_adapters.txt": "stats_adapters_id",
                    },
                },
            },
            "post_trimmed": {
                "{sample}": {
                    "{library}": {
                        "{run}": "stats_fastqc_post",
                    },
                },
            },
        },
    },
}


########################################################################################


class Genome:
    def __init__(self, args, name, filename):
        self.validation_node = ValidateFASTAFilesNode(
            input_file=filename,
            output_file=filename + ".validated",
        )

        # Indexing of FASTA file using 'samtools faidx'
        self.faidx_node = FastaIndexNode(
            infile=filename,
            dependencies=[self.validation_node],
        )

        # Indexing of FASTA file using 'bwa index'
        if args.bwa_algorithm == "mem":
            indexing_class = BWAIndexNode
        elif args.bwa_algorithm == "mem2":
            indexing_class = BWAMem2IndexNode
        else:
            raise RuntimeError(f"unexpected BWA algorithm: {args.bwa_algorithm!r}")

        self.bwa_node = indexing_class(
            input_file=filename,
            dependencies=[self.validation_node],
        )

        # Create sequence dictionary (SAM header)
        self.dict_node = CreateSequenceDictionaryNode(
            in_fasta=filename,
            java_options=args.jre_options,
            dependencies=[self.validation_node],
        )

        # FIXME: Variable number of intervals; disable if 1?
        self.intervals_node = SplitIntervalsNode(
            in_reference=filename,
            out_folder=os.path.join(swap_ext(filename, ".intervals"), "0010"),
            scatter_count=10,
            dependencies=(self.faidx_node, self.dict_node),
        )

        self.name = name
        self.filename = filename
        self.intervals = self.intervals_node.intervals
        self.dependencies = (
            self.faidx_node,
            self.bwa_node,
            self.dict_node,
            self.intervals_node,
        )

    def __str__(self):
        return self.name


def validate_and_index_resources(args, settings):
    # GATK resource files (VCFs); should be tabix indexed
    gatk_resources = set()
    for mode, options in settings["Genotyping"]["VariantRecalibrator"].items():
        if mode in ("SNP", "INDEL"):
            for key, value in options.items():
                if key.startswith("--resource"):
                    gatk_resources.add(value)

    for filename in gatk_resources:
        yield TabixIndexNode(infile=filename, preset="vcf")


def fastqc_sample_runs(args, genome, samples, settings):
    settings = settings["Preprocessing"]
    nodes = []

    # 1. FastQC reports for each input file
    for sample, libraries in samples.items():
        for library, lanes in libraries.items():
            for run, files in lanes.items():
                layout = args.layout.update(
                    sample=sample,
                    library=library,
                    run=run,
                )

                nodes.extend(
                    FastQCNode(
                        in_file=filename,
                        out_folder=layout["stats_fastqc_pre"],
                        options=settings["FastQC"],
                    )
                    for filename in files.values()
                )

                # This mainly serves as a QC step since fastp is too lax
                user_options = {}
                fastp_adapter_options = ("--adapter_sequence", "--adapter_sequence_r2")
                for idx, key in enumerate(fastp_adapter_options, start=1):
                    adapter_sequence = settings["Fastp"].get(key)
                    if adapter_sequence:
                        user_options[f"--adapter{idx}"] = adapter_sequence

                yield IdentifyAdaptersNode(
                    input_file_1=files[1],
                    input_file_2=files[2],
                    output_file=layout["stats_adapters_id"],
                    threads=args.max_threads_fastp,
                    options=user_options,
                )

    # 2. MultiQC report for all files across all samples
    nodes.append(
        MultiQCNode(
            source="fastqc",
            output_prefix=args.layout["stats_fastqc_multiqc_pre"],
            dependencies=nodes,
            options=settings["MultiQC"],
        )
    )

    yield from nodes


def process_fastq_files(args, genome, samples, settings):
    metadata = settings["Metadata"]
    settings = settings["Preprocessing"]

    settings["Fastp"]["--thread"] = args.max_threads_fastp

    nodes = []
    sample = None
    for sample, libraries in samples.items():
        for library, runs in libraries.items():
            for run, files in runs.items():
                layout = args.layout.update(sample=sample, library=library, run=run)

                # Read trimming, quality filtering, and read merging
                fastp_node = FastpNode(
                    in_fq_1=files[1],
                    in_fq_2=files[2],
                    out_fq_1=layout["fastp_paired_1"],
                    out_fq_2=layout["fastp_paired_2"],
                    out_merged=layout["fastp_merged"],
                    out_unpaired=layout["fastp_unpaired"],
                    out_failed=layout["fastp_failed"],
                    out_html=layout["fastp_report_html"],
                    out_json=layout["fastp_report_json"],
                    options=settings["Fastp"],
                )

                # Filtered reads and orphan paired end reads are converted to unmapped
                # BAM alignments so that they can be merged into the final junk BAM
                for key in ("failed", "unpaired"):
                    yield FastqToSamNode(
                        in_fastq=layout["fastp_{}".format(key)],
                        out_bam=layout["fastp_{}_bam".format(key)],
                        options={
                            # Marked as coordinate sorted for faster conversion/merging
                            "--SORT_ORDER": "coordinate",
                            "--READ_GROUP_NAME": run,
                            "--PLATFORM": metadata["Platform"],
                            "--SAMPLE_NAME": sample,
                            "--LIBRARY_NAME": library,
                            "--PLATFORM_UNIT": run,
                            "--DESCRIPTION": "{} reads".format(key.title()),
                        },
                        dependencies=[fastp_node],
                    )

                nodes.append(fastp_node)

                yield fastp_node

    # MultiQC report of reports generated by fastp
    layout = args.layout.update(sample=sample)

    yield MultiQCNode(
        source="fastp",
        output_prefix=layout["stats_fastp_multiqc"],
        options=settings["MultiQC"],
        dependencies=nodes,
    )


def fastqc_trimmed_reads(args, genome, samples, settings):
    settings = settings["Preprocessing"]
    nodes = []

    # 1. FastQC reports for each input file
    for sample, libraries in samples.items():
        for library, lanes in libraries.items():
            for run in lanes:
                layout = args.layout.update(
                    sample=sample,
                    library=library,
                    run=run,
                )

                for read_type in ("fastp_merged", "fastp_paired_1", "fastp_paired_2"):
                    nodes.append(
                        FastQCNode(
                            in_file=layout[read_type],
                            out_folder=layout["stats_fastqc_post"],
                            options=settings["FastQC"],
                        )
                    )

    # 2. MultiQC report for all files across all samples
    nodes.append(
        MultiQCNode(
            source="fastqc",
            output_prefix=args.layout["stats_fastqc_multiqc_post"],
            dependencies=nodes,
            options=settings["MultiQC"],
        )
    )

    return nodes


def map_sample_runs(args, genome, samples, settings):
    metadata = settings["Metadata"]
    bwa_settings = dict(settings["ReadMapping"]["BWAMem"])
    bwa_threads = args.max_threads_bwa

    for sample, libraries in samples.items():
        for library, runs in libraries.items():
            mapped_reads = {
                "merged": [],
                "paired": [],
                "unmapped": [],
            }

            for run in runs:
                layout = args.layout.update(
                    genome=genome,
                    sample=sample,
                    library=library,
                    run=run,
                )

                for read_type in ("fastp_failed_bam", "fastp_unpaired_bam"):
                    mapped_reads["unmapped"].append(layout[read_type])

                filenames = {
                    "paired": (layout["fastp_paired_1"], layout["fastp_paired_2"]),
                    "merged": (layout["fastp_merged"], None),
                }

                for name, (filename_1, filename_2) in filenames.items():
                    layout = layout.update(kind=name)
                    out_bam = layout.get("aln_run_bam", kind=name)
                    mapped_reads[name].append(out_bam)

                    yield BWAAlgorithmNode(
                        reference=genome.filename,
                        input_file_1=filename_1,
                        input_file_2=filename_2,
                        output_file=out_bam,
                        # FIXME: Remove constructor option and only use -t instead?
                        algorithm=args.bwa_algorithm,
                        alt_aware=True,
                        alt_optimize=True,
                        threads=bwa_threads,
                        mapping_options=bwa_settings,
                        # Options passed to 'paleomix cleanup'
                        cleanup_options={
                            # Add mate-score tag required by `samtools markdup`
                            "--add-mate-score": None,
                            "--rg-id": run,
                            "--rg": [
                                "PL:" + metadata["Platform"],
                                "SM:" + sample,
                                "LB:" + library,
                                "PU:" + run,
                                "DS:" + name.title() + " reads",
                            ],
                        },
                    )

            libraries[library] = mapped_reads


def filter_pcr_duplicates(args, genome, samples, settings):
    mode = settings["ReadMapping"]["PCRDuplicates"]["mode"]
    if mode == "skip":
        for sample, libraries in samples.items():
            for library, read_types in libraries.items():
                bam_files = []
                for library_files in read_types.values():
                    bam_files.extend(library_files)

                libraries[library] = bam_files

        return
    elif mode not in ("mark", "filter"):
        raise RuntimeError(f"unknown PCRDuplicates mode {mode!r}")

    markdup_options = {}
    if mode == "filter":
        markdup_options["-r"] = None
    if args.max_threads_samtools > 1:
        markdup_options["--threads"] = args.max_threads_samtools

    for sample, libraries in samples.items():
        for library, read_types in libraries.items():
            layout = args.layout.update(genome=genome, sample=sample, library=library)

            yield MarkDupNode(
                in_bams=read_types["paired"],
                out_bam=layout["aln_rmdup_paired_bam"],
                out_stats=layout["aln_rmdup_paired_metrics"],
                options=markdup_options,
            )

            yield FilterCollapsedBAMNode(
                input_bams=read_types["merged"],
                output_bam=layout["aln_rmdup_merged_bam"],
                keep_dupes=mode == "mark",
            )

            libraries[library] = [
                layout["aln_rmdup_paired_bam"],
                layout["aln_rmdup_merged_bam"],
            ]

            libraries[library].extend(read_types["unmapped"])


def merge_samples_alignments(args, genome, samples, settings):
    for sample, libraries in samples.items():
        input_libraries = []
        for library in libraries.values():
            input_libraries.extend(library)

        layout = args.layout.update(genome=genome, sample=sample)

        # Split BAM into file containing proper alignment and BAM containing junk
        split = FinalizeBAMNode(
            in_bams=input_libraries,
            out_passed=layout["aln_split_passed_bam"],
            out_failed=layout["aln_split_failed_bam"],
            out_json=layout["aln_split_statistics"],
            threads=args.max_threads_samtools,
        )

        samples[sample] = BAMIndexNode(
            infile=layout["aln_split_passed_bam"],
            dependencies=[split],
            options={
                "-@": args.max_threads_samtools,
            },
        )

        yield samples[sample]


def recalibrate_nucleotides(args, genome, samples, settings):
    settings = settings["ReadMapping"]
    # FIXME: Nicer way to specify known sites. Maybe just remove known_sites param?
    recalibrator_options = dict(settings["BaseRecalibrator"])
    recalibrator_known_sites = recalibrator_options.pop("--known-sites")

    for sample in samples:
        layout = args.layout.update(genome=genome, sample=sample)

        model = BaseRecalibratorNode(
            in_reference=genome.filename,
            in_known_sites=recalibrator_known_sites,
            in_bam=layout["aln_split_passed_bam"],
            out_table=layout["aln_recal_training_table"],
            out_log=layout["aln_recal_training_log"],
            options=recalibrator_options,
            java_options=args.jre_options,
        )

        yield ApplyBQSRNode(
            in_node=model,
            out_bam=layout["aln_recal_bam"],
            out_log=layout["aln_recal_log"],
            options=settings["ApplyBQSR"],
            java_options=args.jre_options,
            dependencies=[model],
        )


def final_bam_stats(args, genome, samples, settings):
    fastqc_nodes = []
    for sample in samples:
        for method in BAMStatsNode.METHODS:
            layout = args.layout.update(genome=genome, sample=sample, method=method)

            yield BAMStatsNode(
                method=method,
                infile=layout["aln_recal_bam"],
                outfile=layout["bam_stats"],
                options={
                    # Reasonable performance gains from using up to 3-4 threads
                    "--threads": 3,
                },
            )

        layout = args.layout.update(genome=genome, sample=sample)

        # FastQC of proper alignments
        fastqc_nodes.append(
            FastQCNode(
                in_file=layout["aln_recal_bam"],
                out_folder=layout["bam_fastqc_dir"],
            )
        )

        # FastQC of unmapped/filtered reads
        fastqc_nodes.append(
            FastQCNode(
                in_file=layout["aln_split_failed_bam"],
                out_folder=layout["bam_fastqc_dir"],
            )
        )

    yield MultiQCNode(
        source="fastqc",
        output_prefix=args.layout["bam_multiqc_prefix"],
        dependencies=fastqc_nodes,
    )


def haplotype_samples(args, genome, samples, settings):
    settings = settings["Genotyping"]

    sample_gvcfs = {}
    for sample in samples:
        layout = args.layout.update(genome=genome, sample=sample)
        out_vcf = layout["gvcf_per_sample"]

        sample_gvcfs[out_vcf] = HaplotypeCallerNode(
            in_reference=genome.filename,
            in_bam=layout["aln_recal_bam"],
            out_vcf=out_vcf,
            out_log=layout["gvcf_per_sample_log"],
            options=settings["HaplotypeCaller"],
            java_options=args.jre_options,
        )

    layout = args.layout.update(genome=genome)

    if len(sample_gvcfs) == 1 and len(genome.intervals) == 1:
        haplotyping_func = _haplotype_1_sample_1_interval
    elif len(sample_gvcfs) == 1:
        haplotyping_func = _haplotype_1_sample_n_intervals
    elif len(genome.intervals) == 1:
        haplotyping_func = _haplotype_n_samples_1_interval
    else:
        haplotyping_func = _haplotype_n_samples_n_interval

    yield haplotyping_func(
        args=args,
        genome=genome,
        settings=settings,
        layout=layout,
        gvcfs=sample_gvcfs,
    )


def _haplotype_1_sample_1_interval(args, genome, settings, layout, gvcfs):
    # In the simplest case, we only have a single GVCFs. No additional work is
    # needed to obtain a "merged" file.
    (in_gvcf,) = gvcfs

    return GenotypeGVCFs(
        in_reference=genome.filename,
        in_gvcf=in_gvcf,
        out_vcf=layout["vcf_merged"],
        options=settings["GenotypeGVCFs"],
        java_options=args.jre_options,
        dependencies=gvcfs.values(),
    )


def _haplotype_n_samples_1_interval(args, genome, settings, layout, gvcfs):
    task = CombineGVCFsNode(
        in_reference=genome.filename,
        in_variants=gvcfs,
        out_vcf=layout["gvcf_merged"],
        java_options=args.jre_options,
        dependencies=gvcfs.values(),
    )

    return GenotypeGVCFs(
        in_reference=genome.filename,
        in_gvcf=layout["gvcf_merged"],
        out_vcf=layout["vcf_merged"],
        options=settings["GenotypeGVCFs"],
        java_options=args.jre_options,
        dependencies=[task],
    )


def _haplotype_1_sample_n_intervals(args, genome, settings, layout, gvcfs):
    (in_gvcf,) = gvcfs
    tasks = []
    vcfs = []

    for interval in genome.intervals:
        layout = args.layout.update(genome=genome, part=interval["name"])

        options = dict(settings["GenotypeGVCFs"])
        options["--intervals"] = InputFile(interval["filename"])

        task = GenotypeGVCFs(
            in_reference=genome.filename,
            in_gvcf=in_gvcf,
            out_vcf=layout["vcf_merged_part"],
            options=options,
            java_options=args.jre_options,
            dependencies=gvcfs.values(),
        )

        tasks.append(task)
        vcfs.append(layout["vcf_merged_part"])

    task = GatherVcfsNode(
        # Note that in_vcfs must be in genomic order
        in_vcfs=vcfs,
        out_vcf=layout["vcf_merged"],
        java_options=args.jre_options,
        dependencies=tasks,
    )

    return TabixIndexNode(layout["vcf_merged"], dependencies=[task])


def _haplotype_n_samples_n_interval(args, genome, settings, layout, gvcfs):
    tasks = []
    vcfs = []
    for interval in genome.intervals:
        layout = args.layout.update(genome=genome, part=interval["name"])

        task = CombineGVCFsNode(
            in_reference=genome.filename,
            in_variants=gvcfs,
            out_vcf=layout["gvcf_merged_part"],
            options={
                "--intervals": InputFile(interval["filename"]),
                "--ignore-variants-starting-outside-interval": "true",
            },
            java_options=args.jre_options,
            dependencies=gvcfs.values(),
        )

        task = GenotypeGVCFs(
            in_reference=genome.filename,
            in_gvcf=layout["gvcf_merged_part"],
            out_vcf=layout["vcf_merged_part"],
            options=settings["GenotypeGVCFs"],
            java_options=args.jre_options,
            dependencies=[task],
        )

        tasks.append(task)
        vcfs.append(layout["vcf_merged_part"])

    task = GatherVcfsNode(
        # Note that in_vcfs must be in genomic order
        in_vcfs=vcfs,
        out_vcf=layout["vcf_merged"],
        java_options=args.jre_options,
        dependencies=tasks,
    )

    return TabixIndexNode(layout["vcf_merged"], dependencies=[task])


def recalibrate_haplotype(args, genome, samples, settings):
    settings = settings["Genotyping"]
    layout = args.layout.update(genome=genome)

    if settings["VariantRecalibrator"]["Enabled"]:
        # 1a. Build model for SNP recalibration
        snp_node = VariantRecalibratorNode(
            mode="SNP",
            in_reference=genome.filename,
            in_variant=layout["vcf_merged"],
            out_recal=layout["vcf_recal_training_snp_vcf"],
            out_tranches=layout["vcf_recal_training_snp_trances"],
            out_r_plot=layout["vcf_recal_training_snp_r"],
            out_log=layout["vcf_recal_training_snp_vcf_log"],
            options=settings["VariantRecalibrator"]["SNP"],
            java_options=args.jre_options,
        )

        # Custom tranche plot/table
        plot_node = TranchesPlotsNode(
            input_table=layout["vcf_recal_training_snp_trances"],
            output_prefix=layout["vcf_recal_training_snp_trances_extras"],
            dependencies=[snp_node],
        )

        # 1b. Build model for INDEL recalibration
        indel_node = VariantRecalibratorNode(
            mode="INDEL",
            in_reference=genome.filename,
            in_variant=layout["vcf_merged"],
            out_recal=layout["vcf_recal_training_indel_vcf"],
            out_tranches=layout["vcf_recal_training_indel_trances"],
            out_r_plot=layout["vcf_recal_training_indel_r"],
            out_log=layout["vcf_recal_training_indel_vcf_log"],
            options=settings["VariantRecalibrator"]["INDEL"],
            java_options=args.jre_options,
        )

        yield from (snp_node, indel_node, plot_node)

        if settings["ApplyVQSR"]["Enabled"]:
            # 2. Apply SNP recalibration to original BAM
            recal_node = ApplyVQSRNode(
                mode="SNP",
                in_vcf=snp_node.in_variant,
                in_node=snp_node,
                out_vcf=layout["vcf_recal_snp"],
                out_log=layout["vcf_recal_snp_log"],
                options=settings["ApplyVQSR"]["SNP"],
                java_options=args.jre_options,
                dependencies=[snp_node],
            )

            # 3. Apply INDEL recalibration to SNP recalibrated BAM
            yield ApplyVQSRNode(
                mode="INDEL",
                in_vcf=recal_node.out_vcf,
                in_node=indel_node,
                out_vcf=layout["vcf_recal_snp_indel"],
                out_log=layout["vcf_recal_snp_indel_log"],
                options=settings["ApplyVQSR"]["INDEL"],
                java_options=args.jre_options,
                dependencies=[indel_node, recal_node],
            )


def build_pipeline(args, project):
    pipeline = []

    samples = project["Samples"]
    settings = project["Settings"]
    genome = project["Genome"]

    # FIXME: Maybe allow different tools to have different options
    args.jre_options = list(settings["JavaOptions"]) + [
        # Performance profile for longer running tasks
        "-server",
        # Default temporary directory used by individual commands
        "-Djava.io.tmpdir=%(TEMP_DIR)s",
        # Disable graphical components
        "-Djava.awt.headless=true",
        # Single threaded garbage collection to make CPU usage predictable
        "-XX:+UseSerialGC",
    ]

    args.layout = Layout({"{root}": _LAYOUT}, root=args.output)

    # 1. Validate and process genome
    genome = Genome(args, genome["Name"], genome["Path"])
    pipeline.extend(genome.dependencies)

    # 2. Validate and index resource files
    pipeline.extend(validate_and_index_resources(args, settings))

    if not samples:
        logger = logging.getLogger(__name__)
        logger.warning(
            "Project does not contain any samples; genomes will be prepared for later "
            "use, but no mapping/genotyping will be performed"
        )

        return pipeline
    elif args.run_until == "indexing":
        return pipeline

    analytical_steps = {
        # 3. Do quality analysis of input FASTQ files
        "pre-trimming-qc": fastqc_sample_runs,
        # 4. Process FASTQ files to produce mapping-ready reads
        "read-trimming": process_fastq_files,
        # 5. Do quality analysis of trimmed FASTQ files
        "post-trimming-qc": fastqc_trimmed_reads,
        # 6. Map merged runs to genome and genotype samples
        "read-mapping": map_sample_runs,
        # 7. Filter (mark) PCR duplicates for merged and paired reads
        "pcr-duplicate-filtering-0": filter_pcr_duplicates,
        # 8. Merged PCR duplicate filtered libraries to produce an intermediate BAM.
        "pcr-duplicate-filtering": merge_samples_alignments,
        # 9. Recalibrate base qualities using known variable sites
        "base-recalibration": recalibrate_nucleotides,
        # 10. Collect statistics for the final, processed BAM files
        "mapping-statistics": final_bam_stats,
        # 11. Call haplotypes for each sample
        "haplotyping": haplotype_samples,
        # 12. Recalibrate haplotype qualities using known variants
        "haplotype-recalibration": recalibrate_haplotype,
    }

    include_steps = args.run_from == "indexing"
    for analytical_step, func in analytical_steps.items():
        if include_steps or args.run_from == analytical_step:
            pipeline.extend(func(args, genome, samples, settings))
            include_steps = True

        if args.run_until == analytical_step:
            break

    return pipeline
