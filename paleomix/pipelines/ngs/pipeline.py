#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# ruff: noqa: E501
#
from __future__ import annotations

import logging
import os
from typing import Generator, Literal

from paleomix.common.command import InputFile
from paleomix.common.fileutils import swap_ext
from paleomix.common.layout import Layout
from paleomix.node import Node
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
    TranchesPlotsNode,
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
from paleomix.pipelines.ngs.config import PipelineTarget

########################################################################################
# Project layout

# File structure with leaf values representing labels used to refer a given path
_LAYOUT = {
    "alignments": {
        "{sample}.{genome}.bam": "bam_final_passed",
        "{sample}.{genome}.junk.bam": "bam_final_failed",
    },
    "cache": {
        "samples": {
            "{sample}": {
                "alignments": {
                    "{sample}.{genome}.{library}.{run}.{kind}.bam": "bam_run",
                    "{sample}.{genome}.{library}.rmdup.merged.bam": "bam_rmdup_merged",
                    "{sample}.{genome}.{library}.rmdup.paired.bam": "bam_rmdup_paired",
                    "{sample}.{genome}.{library}.rmdup.paired.metrics.txt": "bam_rmdup_paired_metrics",
                    "{sample}.{genome}.good.bam": "bam_split_passed",
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
                "haplotypes": {
                    "{sample}.{genome}.{part}.g.vcf.gz": "gvcf_sample_part",
                },
            }
        },
        "haplotypes": {
            "{genome}.{part}.g.vcf.gz": "gvcf_merged_part",
            "{genome}.{part}.vcf.gz": "vcf_merged_part",
            "{genome}.{part}.snp.vcf.gz": "vcf_recalibrated_SNP_part",
            "{genome}.{part}.indel.vcf.gz": "vcf_recalibrated_INDEL_part",
        },
    },
    "haplotypes": {
        "{sample}.{genome}.g.vcf.gz": "gvcf_final",
    },
    "genotypes": {
        "{genome}.vcf.gz": "vcf_final",
    },
    "statistics": {
        "alignments_multiQC": "bam_multiqc_prefix",
        "alignments": {
            "{sample}.{genome}.recalibration.txt": "bam_recal_training_table",
            "{sample}.{genome}.mapping.json": "bam_split_statistics",
            "{sample}.{genome}.{method}.txt": "bam_stats",
            "fastqc": "bam_fastqc_dir",
        },
        "genotyping": {
            "{genome}.recalibration.snp.r": "vcf_recal_training_SNP_r",
            "{genome}.recalibration.snp.tranches": "vcf_recal_training_SNP_trances",
            "{genome}.recalibration.indel.r": "vcf_recal_training_INDEL_r",
            "{genome}.recalibration.indel.tranches": "vcf_recal_training_INDEL_trances",
            "{genome}.recalibration.indel.vcf.gz": "vcf_recal_training_INDEL_vcf",
            "{genome}.recalibration.snp.vcf.gz": "vcf_recal_training_SNP_vcf",
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


def validate_and_index_resources(settings):
    # GATK resource files (VCFs); should be tabix indexed
    gatk_resources = set()
    for mode, options in settings["Genotyping"]["VariantRecalibrator"].items():
        if mode in ("SNP", "INDEL"):
            for key, value in options.items():
                if key.startswith("--resource"):
                    gatk_resources.add(value)

    for filename in gatk_resources:
        yield TabixIndexNode(infile=filename, preset="vcf")


def fastqc_sample_runs(args, _genome, samples, settings):
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


def process_fastq_files(args, _genome, samples, settings):
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

                if args.target != PipelineTarget.PREPROCESSING:
                    fastp_node.mark_intermediate_files("*.fastq.gz")

                # Filtered reads and orphan paired end reads are converted to unmapped
                # BAM alignments so that they can be merged into the final junk BAM
                for key in ("failed", "unpaired"):
                    to_sam_task = FastqToSamNode(
                        in_fastq=layout[f"fastp_{key}"],
                        out_bam=layout[f"fastp_{key}_bam"],
                        options={
                            # Marked as coordinate sorted for faster conversion/merging
                            "--SORT_ORDER": "coordinate",
                            "--READ_GROUP_NAME": run,
                            "--PLATFORM": metadata["Platform"],
                            "--SAMPLE_NAME": sample,
                            "--LIBRARY_NAME": library,
                            "--PLATFORM_UNIT": run,
                            "--DESCRIPTION": f"{key.title()} reads",
                        },
                        dependencies=[fastp_node],
                    )

                    if args.target != PipelineTarget.PREPROCESSING:
                        to_sam_task.mark_intermediate_files()

                    yield to_sam_task

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


def fastqc_trimmed_reads(args, _genome, samples, settings):
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
                    out_bam = layout.get("bam_run", kind=name)
                    mapped_reads[name].append(out_bam)

                    task = BWAAlgorithmNode(
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

                    task.mark_intermediate_files()
                    yield task

            libraries[library] = mapped_reads


def filter_pcr_duplicates(args, genome, samples, settings):
    mode = settings["ReadMapping"]["PCRDuplicates"]["mode"]
    if mode == "skip":
        for libraries in samples.values():
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

            task = MarkDupNode(
                in_bams=read_types["paired"],
                out_bam=layout["bam_rmdup_paired"],
                out_stats=layout["bam_rmdup_paired_metrics"],
                options=markdup_options,
            )

            task.mark_intermediate_files()
            yield task

            task = FilterCollapsedBAMNode(
                input_bams=read_types["merged"],
                output_bam=layout["bam_rmdup_merged"],
                keep_dupes=mode == "mark",
            )

            task.mark_intermediate_files()
            yield task

            libraries[library] = [
                layout["bam_rmdup_paired"],
                layout["bam_rmdup_merged"],
            ]

            libraries[library].extend(read_types["unmapped"])


def merge_samples_alignments(args, genome, samples, settings):
    bsqr_enabled = (
        settings["ReadMapping"]["BaseRecalibrator"]["Enabled"]
        and settings["ReadMapping"]["ApplyBQSR"]["Enabled"]
    )

    for sample, libraries in samples.items():
        input_libraries = []
        for library in libraries.values():
            input_libraries.extend(library)

        layout = args.layout.update(genome=genome, sample=sample)

        # Write final BAM directly if BSQR is not enabled
        out_passed = (
            layout["bam_split_passed"] if bsqr_enabled else layout["bam_final_passed"]
        )

        # Split BAM into file containing proper alignment and BAM containing junk
        split = FinalizeBAMNode(
            in_bams=input_libraries,
            out_passed=out_passed,
            out_failed=layout["bam_final_failed"],
            out_json=layout["bam_split_statistics"],
            threads=args.max_threads_samtools,
        )

        samples[sample] = BAMIndexNode(
            infile=out_passed,
            dependencies=[split],
            options={
                "-@": args.max_threads_samtools,
            },
        )

        if bsqr_enabled:
            split.mark_intermediate_files(out_passed)
            samples[sample].mark_intermediate_files()

        yield samples[sample]


def recalibrate_nucleotides(args, genome, samples, settings):
    settings = settings["ReadMapping"]

    train_bsqr_options = dict(settings["BaseRecalibrator"])
    train_bsqr_enabled = train_bsqr_options.pop("Enabled")
    # FIXME: Nicer way to specify known sites. Maybe just remove known_sites param?
    train_bsqr_known_sites = train_bsqr_options.pop("--known-sites")

    apply_bsqr_options = dict(settings["ApplyBQSR"])
    apply_bsqr_enabled = apply_bsqr_options.pop("Enabled")

    if train_bsqr_enabled:
        for sample in samples:
            layout = args.layout.update(genome=genome, sample=sample)

            model = BaseRecalibratorNode(
                in_reference=genome.filename,
                in_known_sites=train_bsqr_known_sites,
                in_bam=layout["bam_split_passed"],
                out_table=layout["bam_recal_training_table"],
                options=train_bsqr_options,
                java_options=args.jre_options,
            )

            if apply_bsqr_enabled:
                yield ApplyBQSRNode(
                    in_node=model,
                    out_bam=layout["bam_final_passed"],
                    options=apply_bsqr_options,
                    java_options=args.jre_options,
                    dependencies=[model],
                )
            else:
                yield model


def final_bam_stats(args, genome, samples, _settings):
    fastqc_nodes = []
    for sample in samples:
        for method in ("stats", "idxstats", "flagstats"):
            layout = args.layout.update(genome=genome, sample=sample, method=method)

            yield BAMStatsNode(
                method=method,
                infile=layout["bam_final_passed"],
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
                in_file=layout["bam_final_passed"],
                out_folder=layout["bam_fastqc_dir"],
            )
        )

        # FastQC of unmapped/filtered reads
        fastqc_nodes.append(
            FastQCNode(
                in_file=layout["bam_final_failed"],
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

    for sample in samples:
        layout = args.layout.update(genome=genome, sample=sample)

        gvcfs: list[str] = []
        for interval in genome.intervals:
            out_gvcf = layout.get("gvcf_sample_part", part=interval["name"])
            gvcfs.append(out_gvcf)

            yield HaplotypeCallerNode(
                in_reference=genome.filename,
                in_bam=layout["bam_final_passed"],
                out_vcf=out_gvcf,
                options={
                    "--intervals": InputFile(interval["filename"]),
                    **settings["HaplotypeCaller"],
                },
                java_options=args.jre_options,
            ).mark_intermediate_files()

        # The final gVCFs are used to ease addition of new samples to old runs
        yield GatherVcfsNode(
            # Note that in_vcfs must be in genomic order
            in_vcfs=gvcfs,
            out_vcf=layout["gvcf_final"],
            java_options=args.jre_options,
        )


def genotype_samples(args, genome, samples, settings):
    layout = args.layout.update(genome=genome)
    settings = settings["Genotyping"]

    vcfs: list[str] = []
    for interval in genome.intervals:
        layout = args.layout.update(
            genome=genome,
            part=interval["name"],
        )

        gvcfs = []
        for sample in samples:
            gvcfs.append(layout.get("gvcf_sample_part", sample=sample))

        # Combine partial per-sample gVCFs into into multi-sample paritial gVCF
        yield CombineGVCFsNode(
            in_reference=genome.filename,
            in_variants=gvcfs,
            out_vcf=layout["gvcf_merged_part"],
            java_options=args.jre_options,
        ).mark_intermediate_files()

        # Genotype partial multi-sample gVCF
        yield GenotypeGVCFs(
            in_reference=genome.filename,
            in_gvcf=layout["gvcf_merged_part"],
            out_vcf=layout["vcf_merged_part"],
            options=settings["GenotypeGVCFs"],
            java_options=args.jre_options,
        ).mark_intermediate_files()

        vcfs.append(layout["vcf_merged_part"])


def recalibrate_genotypes(args, genome, _samples, settings):
    settings = settings["Genotyping"]
    if not settings["VariantRecalibrator"]["Enabled"]:
        return

    layout = args.layout.update(genome=genome)

    yield from _recalibrate_vcf_files(
        args=args,
        genome=genome,
        settings=settings,
        layout=layout,
        mode="SNP",
        in_vcf_key="vcf_merged_part",
        out_vcf_key="vcf_recalibrated_SNP_part",
    )

    # INDEL training and recalibration (training is done in parallel, see below)
    yield from _recalibrate_vcf_files(
        args=args,
        genome=genome,
        settings=settings,
        layout=layout,
        mode="INDEL",
        in_vcf_key="vcf_recalibrated_SNP_part",
        out_vcf_key="vcf_recalibrated_INDEL_part",
    )

    if settings["ApplyVQSR"]["Enabled"]:
        vcfs: list[str] = [
            layout.get("vcf_recalibrated_INDEL_part", part=interval["name"])
            for interval in genome.intervals
        ]

        yield GatherVcfsNode(
            # Note that in_vcfs must be in genomic order
            in_vcfs=vcfs,
            out_vcf=layout["vcf_final"],
            java_options=args.jre_options,
        )


def _recalibrate_vcf_files(
    args,
    genome: Genome,
    settings,
    layout: Layout,
    in_vcf_key: str,
    out_vcf_key: str,
    mode: Literal["SNP", "INDEL"],
) -> Generator[Node]:
    # Training is always done on the "raw" VCF to allow SNP/INDEL training in parallel
    in_variants: list[str] = [
        layout.get("vcf_merged_part", part=interval["name"])
        for interval in genome.intervals
    ]

    # 1. Build models for SNP/INDEL recalibration
    model_node = VariantRecalibratorNode(
        mode=mode,
        in_reference=genome.filename,
        in_variants=in_variants,
        out_recal=layout[f"vcf_recal_training_{mode}_vcf"],
        out_tranches=layout[f"vcf_recal_training_{mode}_trances"],
        out_r_plot=layout[f"vcf_recal_training_{mode}_r"],
        options=settings["VariantRecalibrator"][mode],
        java_options=args.jre_options,
    )

    # 2. Custom tranche plot/table for SNPs/indels
    yield TranchesPlotsNode(
        input_table=layout[f"vcf_recal_training_{mode}_trances"],
        output_prefix=layout[f"vcf_recal_training_{mode}_trances"],
        dependencies=[model_node],
    )

    if settings["ApplyVQSR"]["Enabled"]:
        for interval in genome.intervals:
            layout = layout.update(part=interval["name"])

            # 3. Apply SNP/INDEL recalibration to partial BAMs
            # The VCFs are gathered (see above) and can safely be deleted afterwards
            yield ApplyVQSRNode(
                mode=mode,
                in_vcf=layout.get(in_vcf_key, part=interval["name"]),
                in_node=model_node,
                out_vcf=layout.get(out_vcf_key, part=interval["name"]),
                options=settings["ApplyVQSR"][mode],
                java_options=args.jre_options,
            ).mark_intermediate_files()


def build_pipeline(args, project):
    pipeline = []

    samples = project["Samples"]
    settings = project["Settings"]
    genome = project["Genome"]

    # FIXME: Maybe allow different tools to have different options
    args.jre_options = [
        *settings["JavaOptions"],
        # Performance profile for longer running tasks
        "-server",
        # Default temporary directory used by individual commands
        "-Djava.io.tmpdir=%(TEMP_DIR)s",
        # Disable graphical components
        "-Djava.awt.headless=true",
        # Single threaded garbage collection to make CPU usage predictable
        "-XX:+UseSerialGC",
    ]

    # It is difficult to predict exactly how much memory a given task will use so
    # instead the limits are losened and users are encouraged to use tools like
    # `earlyoom` to avoid running out of memory.
    for value in args.jre_options:
        if value.startswith(("-Xmx", "-XX:MaxRAMPercentage=")):
            break
    else:
        args.jre_options.append("-XX:MaxRAMPercentage=95.0")

    args.layout = Layout({"{root}": _LAYOUT}, root=args.output)

    # 1. Validate and process genome
    genome = Genome(args, genome["Name"], genome["Path"])
    pipeline.extend(genome.dependencies)

    # 2. Validate and index resource files
    pipeline.extend(validate_and_index_resources(settings))

    if not samples:
        logger = logging.getLogger(__name__)
        logger.warning(
            "Project does not contain any samples; genomes will be prepared for later "
            "use, but no mapping/genotyping will be performed"
        )

        return pipeline

    analytical_steps = {
        PipelineTarget.PREPROCESSING: (
            # 3. Do quality analysis of input FASTQ files
            fastqc_sample_runs,
            # 4. Process FASTQ files to produce mapping-ready reads
            process_fastq_files,
            # 5. Do quality analysis of trimmed FASTQ files
            fastqc_trimmed_reads,
        ),
        PipelineTarget.ALIGNMENTS: (
            # 6. Map merged runs to genome and genotype samples
            map_sample_runs,
            # 7. Filter (mark) PCR duplicates for merged and paired reads
            filter_pcr_duplicates,
            # 8. Merged PCR duplicate filtered libraries to produce an intermediate BAM.
            merge_samples_alignments,
            # 9. Recalibrate base qualities using known variable sites
            recalibrate_nucleotides,
            # 10. Collect statistics for the final, processed BAM files
            final_bam_stats,
        ),
        PipelineTarget.HAPLOTYPES: (
            # 11. Call haplotypes for each sample
            haplotype_samples,
        ),
        PipelineTarget.GENOTYPES: (
            # 12. Call genotypes for the combined set of samples
            genotype_samples,
            # 13. Recalibrate haplotype qualities using known variants
            recalibrate_genotypes,
        ),
    }

    for target, functions in analytical_steps.items():
        for func in functions:
            pipeline.extend(func(args, genome, samples, settings))

        if target == args.target:
            break

    return pipeline
