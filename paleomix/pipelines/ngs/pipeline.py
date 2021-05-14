import os

from paleomix.common.fileutils import swap_ext

from paleomix.atomiccmd.command2 import InputFile
from paleomix.nodes.bwa import BWAIndexNode, BWAAlgorithmNode
from paleomix.nodes.commands import FilterCollapsedBAMNode, FinalizeBAMNode
from paleomix.nodes.fastp import FastpNode
from paleomix.nodes.fastqc import FastQCNode
from paleomix.nodes.multiqc import MultiQCNode
from paleomix.nodes.samtools import (
    FastaIndexNode,
    BAMIndexNode,
    BAMStatsNode,
    TabixIndexNode,
)
from paleomix.nodes.validation import ValidateFASTAFilesNode

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
    MarkDuplicatesNode,
    SplitIntervalsNode,
    ValidateBAMNode,
    VariantRecalibratorNode,
)
from paleomix.pipelines.ngs.nodes import (
    TranchesPlotsNode,
)


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
            "haplotypes": {
                "{sample}.{genome}.g.vcf.gz": "gvcf_per_sample",
                "{sample}.{genome}.g.vcf.gz.log": "gvcf_per_sample_log",
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
        "reads": {
            "fastp_multiQC": "stats_fastp_multiqc",
            "fastp": {
                "{sample}": {
                    "{library}": {
                        "{run}.html": "fastp_report_html",
                        "{run}_fastp.json": "fastp_report_json",
                    },
                },
            },
            "fastqc_multiQC": "stats_fastqc_multiqc",
            "fastqc": {
                "{sample}": {
                    "{library}": {
                        "{run}": "stats_fastqc",
                    },
                },
            },
        },
    },
}


class Layout:
    def __init__(self, args, **kwargs):
        self.kwargs = {}
        self.kwargs["root"] = args.output

        if "genome" in kwargs:
            kwargs["genome"] = kwargs["genome"].name

        for key, value in kwargs.items():
            if not isinstance(value, str):
                raise ValueError((key, value))

            self.kwargs[key] = value

        if Layout._layout is None:
            Layout._layout = {}
            for key, value in Layout._flatten_layout({"{root}": _LAYOUT}):
                assert key not in Layout._layout, key
                Layout._layout[key] = value

    def __getitem__(self, key):
        return Layout._layout[key].format(**self.kwargs)

    @classmethod
    def _flatten_layout(cls, layout):
        for key, value in layout.items():
            if isinstance(value, str):
                yield value, key
            else:
                for label, path in cls._flatten_layout(value):
                    yield label, os.path.join(key, path)

    _layout = None


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
        self.bwa_node = BWAIndexNode(
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


def fastqc_sample_runs(args, samples, settings):
    nodes = []

    # 1. FastQC reports for each input file
    for sample, libraries in samples.items():
        for library, lanes in libraries.items():
            for run, files in lanes.items():
                layout = Layout(args=args, sample=sample, library=library, run=run)

                nodes.extend(
                    FastQCNode(
                        in_file=filename,
                        out_folder=layout["stats_fastqc"],
                        options=settings["FastQC"],
                    )
                    for filename in files.values()
                )

    # 2. MultiQC report for all files across all samples
    nodes.append(
        MultiQCNode(
            source="fastqc",
            output_prefix=Layout(args)["stats_fastqc_multiqc"],
            dependencies=nodes,
            options=settings["MultiQC"],
        )
    )

    return nodes


def process_fastq_files(args, samples, settings):
    settings_preproc = settings["Preprocessing"]
    settings_meta = settings["Metadata"]

    nodes = []
    for sample, libraries in samples.items():
        for library, runs in libraries.items():
            for run, files in runs.items():
                layout = Layout(args, sample=sample, library=library, run=run)

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
                    options=settings_preproc["Fastp"],
                )

                # Filtered reads and orphan paired end reads are converted to unmapped
                # BAM alignments so that they can be merged into the final junk BAM
                unmapped_reads = {}
                for key in ("failed", "unpaired"):
                    out_bam = layout["fastp_{}_bam".format(key)]
                    unmapped_reads[out_bam] = FastqToSamNode(
                        in_fastq=layout["fastp_{}".format(key)],
                        out_bam=out_bam,
                        options={
                            # Marked as coordinate sorted for faster conversion/merging
                            "--SORT_ORDER": "coordinate",
                            "--READ_GROUP_NAME": run,
                            "--PLATFORM": settings_meta["Platform"],
                            "--SAMPLE_NAME": sample,
                            "--LIBRARY_NAME": library,
                            "--PLATFORM_UNIT": run,
                            "--DESCRIPTION": "{} reads".format(key.title()),
                        },
                        dependencies=[fastp_node],
                    )

                nodes.append(fastp_node)
                runs[run] = {
                    "node": fastp_node,
                    "out_merged_1": layout["fastp_merged"],
                    "out_merged_2": None,
                    "out_paired_1": layout["fastp_paired_1"],
                    "out_paired_2": layout["fastp_paired_2"],
                    "unmapped": unmapped_reads,
                }

    # MultiQC report of reports generated by fastp
    layout = Layout(args, sample=sample)

    yield MultiQCNode(
        source="fastp",
        output_prefix=layout["stats_fastp_multiqc"],
        options=settings_preproc["MultiQC"],
        dependencies=nodes,
    )


def map_sample_runs(args, genome, samples, settings):
    # FIXME: Should be handled by Node itself
    bwa_settings = dict(settings["BWAMem"])
    bwa_threads = bwa_settings.pop("-t", 4)

    for sample, libraries in samples.items():
        for library, runs in libraries.items():
            mapped_reads = {
                "merged": {},
                "paired": {},
                "unmapped": {},
            }

            for run, files in runs.items():
                # Unmapped reads
                mapped_reads["unmapped"].update(files["unmapped"])

                for name in ("merged", "paired"):
                    layout = Layout(
                        args=args,
                        genome=genome,
                        sample=sample,
                        library=library,
                        run=run,
                        kind=name,
                    )

                    out_bam = layout["aln_run_bam"]
                    bwa = BWAAlgorithmNode(
                        reference=genome.filename,
                        input_file_1=files["out_{}_1".format(name)],
                        input_file_2=files["out_{}_2".format(name)],
                        output_file=out_bam,
                        # FIXME: Remove constructor option and only use -t instead?
                        threads=bwa_threads,
                        mapping_options=bwa_settings,
                        # Options passed to 'paleomix cleanup'
                        cleanup_options={
                            "--rg-id": run,
                            "--rg": [
                                # FIXME: PL should not be hardcoded
                                "PL:ILLUMINA",
                                "SM:" + sample,
                                "LB:" + library,
                                "PU:" + run,
                                "DS:" + name.title() + " reads",
                            ],
                        },
                        dependencies=[files["node"]],
                    )

                    mapped_reads[name][out_bam] = ValidateBAMNode(
                        in_bam=out_bam,
                        java_options=args.jre_options,
                        out_log=layout["aln_run_validation_log"],
                        dependencies=[bwa],
                    )

            libraries[library] = mapped_reads

    return ()


def filter_pcr_duplicates(args, genome, samples, settings):
    for sample, libraries in samples.items():
        for library, read_types in libraries.items():
            layout = Layout(args, genome=genome, sample=sample, library=library)

            paired = MarkDuplicatesNode(
                in_bams=read_types["paired"].keys(),
                out_bam=layout["aln_rmdup_paired_bam"],
                out_metrics=layout["aln_rmdup_paired_metrics"],
                java_options=args.jre_options,
                dependencies=read_types["paired"].values(),
            )

            merged = FilterCollapsedBAMNode(
                input_bams=read_types["merged"].keys(),
                output_bam=layout["aln_rmdup_merged_bam"],
                keep_dupes=True,
                dependencies=read_types["merged"].values(),
            )

            libraries[library] = {
                layout["aln_rmdup_paired_bam"]: paired,
                layout["aln_rmdup_merged_bam"]: merged,
            }

            libraries[library].update(read_types["unmapped"])

    return ()


def merge_samples_alignments(args, genome, samples, settings):
    for sample, libraries in samples.items():
        input_libraries = {}
        for library in libraries.values():
            input_libraries.update(library)

        layout = Layout(args, genome=genome, sample=sample)

        # Split BAM into file containing proper alignment and BAM containing junk
        split = FinalizeBAMNode(
            in_bams=input_libraries.keys(),
            out_passed=layout["aln_split_passed_bam"],
            out_failed=layout["aln_split_failed_bam"],
            out_json=layout["aln_split_statistics"],
            options={
                # Reasonable performance gains from using up to 3-4 threads
                "--threads": 3,
            },
            dependencies=input_libraries.values(),
        )

        indexed = BAMIndexNode(
            infile=layout["aln_split_passed_bam"],
            dependencies=[split],
            options={
                # Reasonable performance gains from using up to 3-4 threads
                "-@": 3,
            },
        )

        samples[sample] = ValidateBAMNode(
            in_bam=layout["aln_split_passed_bam"],
            in_index=layout["aln_split_passed_bam"] + ".bai",
            java_options=args.jre_options,
            out_log=layout["aln_split_good_validation_log"],
            dependencies=[indexed],
        )

    return ()


def recalibrate_nucleotides(args, genome, samples, settings):
    # FIXME: Nicer way to specify known sites. Maybe just remove known_sites param?
    recalibrator_options = dict(settings["BaseRecalibrator"])
    recalibrator_known_sites = recalibrator_options.pop("--known-sites")

    for sample, node in samples.items():
        layout = Layout(args, genome=genome, sample=sample)

        model = BaseRecalibratorNode(
            in_reference=genome.filename,
            in_known_sites=recalibrator_known_sites,
            in_bam=node.out_bam,
            out_table=layout["aln_recal_training_table"],
            out_log=layout["aln_recal_training_log"],
            options=recalibrator_options,
            java_options=args.jre_options,
            dependencies=[node],
        )

        recalibrator = ApplyBQSRNode(
            in_node=model,
            out_bam=layout["aln_recal_bam"],
            out_log=layout["aln_recal_log"],
            options=settings["ApplyBQSR"],
            java_options=args.jre_options,
            dependencies=[model],
        )

        samples[sample] = recalibrator
        yield recalibrator


def final_bam_stats(args, genome, samples, settings):
    fastqc_nodes = []
    for sample, node in samples.items():
        for method in BAMStatsNode.METHODS:
            layout = Layout(args, genome=genome, sample=sample, method=method)

            yield BAMStatsNode(
                method=method,
                infile=node.out_bam,
                outfile=layout["bam_stats"],
                options={
                    # Reasonable performance gains from using up to 3-4 threads
                    "--threads": 3,
                },
                dependencies=[node],
            )

        layout = Layout(args, genome=genome, sample=sample)

        # FastQC of proper alignments
        fastqc_nodes.append(
            FastQCNode(
                in_file=node.out_bam,
                out_folder=layout["bam_fastqc_dir"],
                dependencies=[node],
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
        output_prefix=layout["bam_multiqc_prefix"],
        dependencies=fastqc_nodes,
    )


def haplotype_samples(args, genome, samples, settings):
    sample_gvcfs = {}
    for sample, node in samples.items():
        layout = Layout(args, genome=genome, sample=sample)
        out_vcf = layout["gvcf_per_sample"]

        sample_gvcfs[out_vcf] = HaplotypeCallerNode(
            in_reference=genome.filename,
            in_bam=node.out_bam,
            out_vcf=out_vcf,
            out_log=layout["gvcf_per_sample_log"],
            options=settings["HaplotypeCaller"],
            java_options=args.jre_options,
            dependencies=[node],
        )

    layout = Layout(args, genome=genome)

    if len(sample_gvcfs) == 1:
        # In the simplest case, we only have a single GVCFs. No additional work is
        # needed to obtain a "merged" file.
        (combined_gvcf,) = sample_gvcfs
        (combined_node,) = sample_gvcfs.values()
    elif len(genome.intervals) > 1:
        part_nodes = []
        for interval in genome.intervals:
            layout = Layout(args, genome=genome, part=interval["name"])

            part_nodes.append(
                CombineGVCFsNode(
                    in_reference=genome.filename,
                    in_variants=sample_gvcfs,
                    out_vcf=layout["gvcf_merged_part"],
                    options={
                        "--intervals": InputFile(interval["filename"]),
                        "--ignore-variants-starting-outside-interval": "true",
                    },
                    java_options=args.jre_options,
                    dependencies=sample_gvcfs.values(),
                )
            )

        # Combine scattered GVCFs into a single
        combined_gvcf = layout["gvcf_merged"]
        gather_node = GatherVcfsNode(
            # Note that in_vcfs must be in genomic order
            in_vcfs=[node.out_vcf for node in part_nodes],
            out_vcf=combined_gvcf,
            java_options=args.jre_options,
            dependencies=part_nodes,
        )

        combined_node = TabixIndexNode(infile=combined_gvcf, dependencies=[gather_node])

        yield from part_nodes
    else:
        combined_gvcf = layout["gvcf_merged"]
        combined_node = CombineGVCFsNode(
            in_reference=genome.filename,
            in_variants=sample_gvcfs,
            out_vcf=combined_gvcf,
            java_options=args.jre_options,
            dependencies=sample_gvcfs.values(),
        )

    yield GenotypeGVCFs(
        in_reference=genome.filename,
        in_gvcf=combined_gvcf,
        out_vcf=layout["vcf_merged"],
        options=settings["GenotypeGVCFs"],
        java_options=args.jre_options,
        dependencies=[combined_node],
    )


def recalibrate_haplotype(args, genome, samples, settings):
    layout = Layout(args, genome=genome)

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
                in_node=snp_node,
                out_vcf=layout["vcf_recal_snp"],
                out_log=layout["vcf_recal_snp_log"],
                options=settings["ApplyVQSR"]["SNP"],
                java_options=args.jre_options,
                dependencies=[snp_node],
            )

            # 3. Apply INDEL recalibration to original SNP recalibrated BAM
            yield ApplyVQSRNode(
                mode="INDEL",
                in_node=indel_node,
                out_vcf=layout["vcf_recal_snp_indel"],
                out_log=layout["vcf_recal_snp_indel_log"],
                options=settings["ApplyVQSR"]["INDEL"],
                java_options=args.jre_options,
                dependencies=[indel_node, recal_node],
            )


def build_pipeline(args, project):
    pipeline = []

    def _add(func, *args_, **kwargs):
        pipeline.extend(func(args, *args_, **kwargs))

    samples = project["Samples"]
    settings = project["Settings"]
    genome = project["Genome"]

    # FIXME: Maybe allow different tools to have different options
    args.jre_options = [
        "-server",
        "-Djava.io.tmpdir={}".format(args.temp_root),
        "-Djava.awt.headless=true",
    ]

    args.jre_options.extend(settings["JavaOptions"])

    # 1. Validate and process genome
    genome = Genome(args, genome["Name"], genome["Path"])
    pipeline.extend(genome.dependencies)

    # 2. Validate and index resource files
    _add(validate_and_index_resources, settings)

    # 3. Do quality analysis of input FASTQ files
    _add(fastqc_sample_runs, samples, settings["Preprocessing"])

    # 4. Process FASTQ files to produce mapping-ready reads
    _add(process_fastq_files, samples, settings)

    # 5. Map merged runs to genome and genotype samples
    _add(map_sample_runs, genome, samples, settings["ReadMapping"])

    # 6. Filter (mark) PCR duplicates for merged and paired reads
    _add(filter_pcr_duplicates, genome, samples, settings["ReadMapping"])

    # 7. Merged PCR duplicate filtered libraries to produce an intermediate BAM.
    _add(merge_samples_alignments, genome, samples, settings["ReadMapping"])

    # 9. Recalibrate base qualities using known variable sites
    _add(recalibrate_nucleotides, genome, samples, settings["ReadMapping"])

    # 9. Collect statistics for the final, processed BAM files
    _add(final_bam_stats, genome, samples, settings["ReadMapping"])

    # 10. Call haplotypes for each sample
    _add(haplotype_samples, genome, samples, settings["Genotyping"])

    # 11. Recalibrate haplotype qualities using known variants
    _add(recalibrate_haplotype, genome, samples, settings["Genotyping"])

    return pipeline
