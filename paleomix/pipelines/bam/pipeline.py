#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import logging
import os
from collections import defaultdict

import paleomix
import paleomix.common.logging
import paleomix.node
from paleomix.common.fileutils import swap_ext
from paleomix.common.layout import Layout
from paleomix.common.yaml import YAMLError
from paleomix.nodes.adapterremoval import PE_AdapterRemovalNode, SE_AdapterRemovalNode
from paleomix.nodes.bowtie2 import Bowtie2IndexNode, Bowtie2Node
from paleomix.nodes.bwa import (
    BWAAlgorithmNode,
    BWABacktrack,
    BWAIndexNode,
    BWASampe,
    BWASamse,
)
from paleomix.nodes.commands import FilterCollapsedBAMNode
from paleomix.nodes.mapdamage import (
    MapDamageModelNode,
    MapDamagePlotNode,
    MapDamageRescaleNode,
)
from paleomix.nodes.samtools import (
    BAMIndexNode,
    BAMMergeNode,
    FastaIndexNode,
    MarkDupNode,
)
from paleomix.nodes.validation import ValidateFASTAFilesNode, ValidateFASTQFilesNode
from paleomix.pipeline import Pypeline
from paleomix.pipelines.bam.makefile import MakefileError, read_makefiles

LAYOUT = {
    "{sample}.cache": {
        "{genome}.validated": "final_validated",
        "{genome}": {
            "{library}": {
                "{lane}": {
                    "{shortname}": {
                        "{read_type}.bam": "initial_bam",
                        "{read_type}.stats": "initial_stats",
                    }
                }
            },
            # Libraries processed using `samtools markdup` or `paleomix rmdup_collapsed`
            "{library}.rmdup.{method}.bam": "deduplicated_bam",
            "{library}.rmdup.{method}.statistics": "deduplicated_stats",
            # Libraries where quality scores have been rescaled using mapDamage
            "{library}.rescaled.bam": "rescaled_bam",
        },
        "reads": {
            "{library}": {
                "{lane}": {
                    "{shortname}": {
                        "reads": "reads_prefix",
                        "pretrimmed.json": "reads_statistics",
                    },
                },
            },
        },
    },
    "{sample}.{genome}.bam": "final_bam",
    "{sample}.{genome}.bam.bai": "final_bai",
    "{sample}.{genome}.mapDamage": {
        "{library}": "mapdamage_folder",
    },
}

########################################################################################


def index_genomes(log, makefiles):
    tasks = {}
    any_errors = False
    for makefile in makefiles:
        genomes = makefile["Genomes"]
        if not genomes:
            log.error("No genomes specified in %r", makefile["Filename"])
            any_errors = True

        for genome in genomes.values():
            path = genome["Path"]
            # Multiple genomes may use the path from different
            abspath = os.path.abspath(path)

            if abspath not in tasks:
                # Basic validation is preformed since downstream tools may produce
                # inconsistent/unexpected results when run on malformed files, such
                # as `samtools faidx` ignoring any sequence with a previously used name.
                validation = ValidateFASTAFilesNode(
                    input_file=path,
                    output_file=path + ".validated",
                )

                # Indexing of FASTA file using 'samtools faidx'
                fai_indexing = FastaIndexNode(path, dependencies=[validation])

                # Indexing of FASTA file using 'bwa index'
                bwa_indexing = BWAIndexNode(path, dependencies=[validation])

                # Indexing of FASTA file using 'bowtie2-build'
                bowtie2_indexing = Bowtie2IndexNode(path, dependencies=[validation])

                tasks[abspath] = {
                    "BWA": [fai_indexing, bwa_indexing],
                    "Bowtie2": [fai_indexing, bowtie2_indexing],
                }

            genome["Tasks"] = tasks[abspath]

    return not any_errors


########################################################################################


def process_fastq_reads(args, layout, record):
    lane_type = record["Type"]
    if lane_type == "Untrimmed":
        trimmed_reads = _process_untrimmed_reads(layout, record, args)
    else:
        trimmed_reads = _process_pretrimmed_reads(layout, record)

    records = []
    for read_type, files, task in trimmed_reads:
        # TODO: Warn if lane is completely excluded
        if not record["Options"]["ExcludeReads"][read_type]:
            records.append(
                {
                    "Task": task,
                    "Path": files,
                    "Type": read_type,
                    "Shortname": record["Shortname"],
                    "Options": record["Options"],
                }
            )

    return records


def _process_pretrimmed_reads(layout, record):
    read_type = record["Type"]
    input_files = [filename for filename in record["Path"] if filename is not None]

    yield read_type, record["Path"], ValidateFASTQFilesNode(
        input_files=input_files,
        output_file=layout.get(
            "reads_statistics",
            shortname=record["Shortname"],
        ),
        offset=record["Options"]["QualityOffset"],
        collapsed=("Collapsed" in read_type),
    )


def _process_untrimmed_reads(layout, record, args):
    quality_offset = record["Options"]["QualityOffset"]
    options = dict(record["Options"]["AdapterRemoval"])

    if quality_offset != 33:
        options["--qualitybase"] = quality_offset
        # Quality scores of trimmed reads is normalized to Phred+33
        options["--qualitybase-output"] = 33

    output_prefix = layout.get(
        "reads_prefix",
        shortname=record["Shortname"],
    )

    file_1, file_2 = record["Path"]

    if file_2 is None:
        task = SE_AdapterRemovalNode(
            input_file=file_1,
            output_prefix=output_prefix,
            threads=args.adapterremoval_max_threads,
            options=options,
        )

        yield ("Single", (task.out_truncated, None), task)
    else:
        task = PE_AdapterRemovalNode(
            input_file_1=file_1,
            input_file_2=file_2,
            output_prefix=output_prefix,
            threads=args.adapterremoval_max_threads,
            options=options,
        )

        yield ("Singleton", (task.out_singleton, None), task)

        out_paired = (task.out_paired.format(Pair=1), task.out_paired.format(Pair=2))
        yield ("Paired", out_paired, task)

        if task.out_merged:
            yield ("Collapsed", (task.out_merged, None), task)

        if task.out_merged_truncated:
            yield ("CollapsedTruncated", (task.out_merged_truncated, None), task)


########################################################################################


def map_fastq_reads(args, layout, genome, record):
    layout = layout.update(
        shortname=record["Shortname"],
        read_type=record["Type"].lower(),
    )

    options = record["Options"]
    aligner = options["Aligners"]["Program"]
    output_file = layout["initial_bam"]

    # Common mapping parameters
    parameters = {
        "input_file_1": record["Path"][0],
        "input_file_2": record["Path"][1],
        "output_file": output_file,
        "reference": genome["Path"],
        "mapping_options": record["Options"]["Aligners"][aligner],
        "cleanup_options": _cleanup_options(record, layout),
        "dependencies": [record["Task"]] + genome["Tasks"][aligner],
    }

    if aligner == "BWA":
        algorithm = options["Aligners"][aligner]["Algorithm"].lower()
        parameters["threads"] = args.bwa_max_threads

        if algorithm == "backtrack":
            mapping_task_func = _build_bwa_backtrack_task
        elif algorithm in ("mem", "bwasw"):
            parameters["algorithm"] = algorithm

            mapping_task_func = BWAAlgorithmNode
        else:
            raise NotImplementedError(f"BWA {algorithm} not implemented!")
    elif aligner == "Bowtie2":
        parameters["threads"] = args.bowtie2_max_threads

        mapping_task_func = Bowtie2Node
    else:
        raise NotImplementedError(f"Aligner {aligner!r} not supported!")

    return {
        "Type": record["Type"],
        "Path": output_file,
        "Options": record["Options"],
        "Task": mapping_task_func(**parameters),
    }


def _build_bwa_backtrack_se_task(
    *,
    input_file,
    output_file,
    reference,
    threads,
    dependencies,
    mapping_options,
    cleanup_options,
):
    output_file_sai = swap_ext(output_file, ".sai")

    sai_task = BWABacktrack(
        input_file=input_file,
        output_file=output_file_sai,
        threads=threads,
        reference=reference,
        mapping_options=mapping_options,
        dependencies=dependencies,
    )

    return BWASamse(
        input_file_fq=input_file,
        input_file_sai=output_file_sai,
        output_file=output_file,
        reference=reference,
        threads=max(2, threads // 2),
        cleanup_options=cleanup_options,
        dependencies=sai_task,
    )


def _build_bwa_backtrack_pe_task(
    *,
    input_file_1,
    input_file_2,
    output_file,
    reference,
    threads,
    dependencies,
    mapping_options,
    cleanup_options,
):
    backtrack_options = {
        "threads": threads,
        "reference": reference,
        "mapping_options": mapping_options,
        "dependencies": dependencies,
    }

    output_sai_1 = swap_ext(output_file, "%i.sai" % (1,))
    output_sai_2 = swap_ext(output_file, "%i.sai" % (2,))

    task_sai_1 = BWABacktrack(input_file_1, output_sai_1, **backtrack_options)
    task_sai_2 = BWABacktrack(input_file_2, output_sai_2, **backtrack_options)

    return BWASampe(
        input_file_sai_1=output_sai_1,
        input_file_sai_2=output_sai_2,
        input_file_fq_1=input_file_1,
        input_file_fq_2=input_file_2,
        output_file=output_file,
        reference=reference,
        threads=max(2, threads // 2),
        cleanup_options=cleanup_options,
        dependencies=(task_sai_1, task_sai_2),
    )


def _build_bwa_backtrack_task(input_file_1, input_file_2, mapping_options, **kwargs):
    if not mapping_options["UseSeed"]:
        mapping_options = dict(mapping_options)
        mapping_options["-l"] = 2 ** 16 - 1

    if input_file_2 is None:
        return _build_bwa_backtrack_se_task(
            input_file=input_file_1,
            mapping_options=mapping_options,
            **kwargs,
        )
    else:
        return _build_bwa_backtrack_pe_task(
            input_file_1=input_file_1,
            input_file_2=input_file_2,
            mapping_options=mapping_options,
            **kwargs,
        )


def _cleanup_options(record, layout):
    aligner = record["Options"]["Aligners"]["Program"]
    aligner_options = record["Options"]["Aligners"][aligner]
    platform = record["Options"]["Platform"].upper()

    sample = layout.get_field("sample")
    library = layout.get_field("library")
    barcode = layout.get_field("lane")

    options = {
        # Add mate-score tag required by `samtools markdup`
        "--add-mate-score": None,
        "--rg-id": library,
        "--rg": [
            f"SM:{sample}",
            f"LB:{library}",
            f"PU:{barcode}",
            f"PL:{platform}",
        ],
        "-q": aligner_options["MinQuality"],
    }

    if aligner_options["FilterUnmappedReads"]:
        options["-F"] = "0x4"

    return options


########################################################################################


def filter_pcr_duplicates(args, layout, records):
    # The PCRDuplicates feature cannot be set below the library level,
    # checking the first record is sufficient
    pcr_filtering = records[0]["Options"]["Features"]["PCRDuplicates"]
    if not pcr_filtering:
        return records

    options = None
    tasks_by_read_type = defaultdict(dict)
    for record in records:
        filepath = record["Path"]
        options = record["Options"]

        if record["Type"] in ("Collapsed", "CollapsedTruncated"):
            tasks_by_read_type["Merged"][filepath] = record["Task"]
        elif record["Type"] in ("Single", "Paired", "Singleton"):
            tasks_by_read_type["Unmerged"][filepath] = record["Task"]
        else:
            raise NotImplementedError(record)

    return _filter_pcr_duplicates_by_type(
        args=args,
        layout=layout,
        options=options,
        tasks_by_read_type=tasks_by_read_type,
        strategy=pcr_filtering,
    )


def _filter_pcr_duplicates_by_type(args, layout, options, tasks_by_read_type, strategy):
    keep_duplicates = isinstance(strategy, str) and strategy.lower() == "mark"
    markdup_options = {"--threads": args.samtools_max_threads}
    if not keep_duplicates:
        markdup_options["-r"] = None

    records = []
    for key, filenames_and_tasks in tasks_by_read_type.items():
        layout = layout.update(method=key.lower())
        out_bam = layout["deduplicated_bam"]
        out_statistics = layout["deduplicated_stats"]

        if key == "Merged":
            task = FilterCollapsedBAMNode(
                input_bams=list(filenames_and_tasks),
                output_bam=out_bam,
                keep_dupes=keep_duplicates,
                dependencies=filenames_and_tasks.values(),
            )
        elif key == "Unmerged":
            task = MarkDupNode(
                in_bams=list(filenames_and_tasks),
                out_bam=out_bam,
                out_stats=out_statistics,
                options=markdup_options,
                dependencies=filenames_and_tasks.values(),
            )
        else:
            raise RuntimeError("unexpected read type {!r}".format(key))

        records.append(
            {
                "Type": key,
                "Path": out_bam,
                "Task": task,
                "Options": options,
            }
        )

    return records


########################################################################################


def run_mapdamage(layout: Layout, genome, records):
    options = records[0]["Options"]
    run_type = options["Features"]["mapDamage"]

    extra_task = None
    if run_type in ("rescale", "model", "plot", True):
        # Basic run of mapDamage, only generates plots / tables
        extra_task = MapDamagePlotNode(
            reference=genome["Path"],
            input_files=[record["Path"] for record in records],
            output_directory=layout["mapdamage_folder"],
            title="mapDamage plot for library %r" % (layout.get_field("library"),),
            options=options["mapDamage"],
            dependencies=[record["Task"] for record in records],
        )

    if run_type in ("rescale", "model"):
        # Builds model of post-mortem DNA damage
        assert extra_task is not None
        extra_task = MapDamageModelNode(
            reference=genome["Path"],
            directory=layout["mapdamage_folder"],
            options=options,
            dependencies=(extra_task,),
        )

    if run_type in ("rescale",):
        # Rescales BAM quality scores using model built above
        assert extra_task is not None
        task = MapDamageRescaleNode(
            reference=genome["Path"],
            input_files=[record["Path"] for record in records],
            output_file=layout["rescaled_bam"],
            directory=layout["mapdamage_folder"],
            options=options["mapDamage"],
            dependencies=(extra_task,),
        )

        extra_task = None
        records = [
            {
                "Type": "Rescaled",
                "Path": layout["rescaled_bam"],
                "Task": task,
                "Options": options,
            }
        ]

    return records, extra_task


########################################################################################


def merge_libraries(args, layout, genome, records):
    # FIXME: Do a file-copy if there is only one input file
    task = BAMMergeNode(
        in_files=[record["Path"] for record in records],
        out_file=layout["final_bam"],
        options={
            "--threads": args.samtools_max_threads,
        },
        dependencies=[record["Task"] for record in records],
    )

    task = BAMIndexNode(
        infile=layout["final_bam"],
        index_format=genome["IndexFormat"],
        options={
            "-@": args.samtools_max_threads,
        },
        dependencies=[task],
    )

    return {
        "Genome": genome,
        "Path": layout["final_bam"],
        "Task": task,
        "Options": records[0]["Options"],
    }


########################################################################################


def build_pipeline_trimming(args, makefile):
    for sample, libraries in makefile["Samples"].items():
        for library, lanes in libraries.items():
            for barcode, records in lanes.items():
                for record in records:
                    layout = args.layout.update(
                        sample=sample,
                        library=library,
                        lane=barcode,
                    )

                    for trimmed_reads in process_fastq_reads(args, layout, record):
                        yield trimmed_reads["Task"]


def build_pipeline_full(args, makefile):
    layout = args.layout
    for sample, library_records in makefile["Samples"].items():
        # Trimmed reads are instantiated as needed, ensuring that sorting tasks by ID
        # results in trimming (IO heavy) and mapping (CPU heavy) tasks being interleaved
        trimmed_reads_cache = {}

        for genome in makefile["Genomes"].values():
            libraries = []

            for library, lane_records in library_records.items():
                lanes = []

                for barcode, records in lane_records.items():
                    layout = args.layout.update(
                        sample=sample,
                        genome=genome["Name"],
                        library=library,
                        lane=barcode,
                    )

                    for record in records:
                        key = id(record)
                        trimmed_reads = trimmed_reads_cache.get(key)
                        if trimmed_reads is None:
                            # Trim untrimmed reads and validate already trimmed reads
                            trimmed_reads = process_fastq_reads(args, layout, record)
                            trimmed_reads_cache[key] = trimmed_reads

                        for record in trimmed_reads:
                            # Map trimmed reads to the target genome
                            lanes.append(map_fastq_reads(args, layout, genome, record))

                # Exclusion of read types may result in the elimination of libraries
                if lanes:
                    # Optionally filter PCR duplicates
                    library = filter_pcr_duplicates(args, layout, lanes)
                    # Optionally run mapDamage on library
                    library, extra_task = run_mapdamage(layout, genome, library)
                    if extra_task is not None:
                        yield extra_task

                    libraries.extend(library)

            # Exclusion of read types may lead to the elimination of all libraries
            if libraries:
                record = merge_libraries(args, layout, genome, libraries)

                yield record["Task"]

        # TODO: Per sample statistics


def run(config, pipeline_variant):
    paleomix.common.logging.initialize(
        log_level=config.log_level,
        log_file=config.log_file,
        auto_log_file=os.path.join(config.temp_root, "bam_pipeline"),
    )

    logger = logging.getLogger(__name__)
    if pipeline_variant not in ("bam", "trim"):
        logger.critical("Unexpected BAM pipeline variant %r", pipeline_variant)
        return 1

    if not os.path.exists(config.temp_root):
        try:
            os.makedirs(config.temp_root)
        except OSError as error:
            logger.error("Could not create temp root: %s", error)
            return 1

    if not os.access(config.temp_root, os.R_OK | os.W_OK | os.X_OK):
        logger.error("Insufficient permissions for temp root: %r", config.temp_root)
        return 1

    try:
        makefiles = read_makefiles(config.makefiles)
    except (MakefileError, YAMLError, IOError) as error:
        logger.error("Error reading makefiles: %s", error)
        return 1

    pipeline_func = build_pipeline_trimming
    if pipeline_variant != "trim":
        # Genomes are processed first so that these tasks are started before reads are
        # trimmed. This allows mapping to be started as early as possible.
        if not index_genomes(logger, makefiles):
            return 1

        pipeline_func = build_pipeline_full

    config.layout = Layout({"{root}": LAYOUT}, root=config.destination)

    nodes = []
    for makefile in makefiles:
        logger.info("Building BAM pipeline for %r", makefile["Filename"])
        try:
            nodes.extend(pipeline_func(config, makefile))
        except paleomix.node.NodeError as error:
            logger.error(
                "Error while building pipeline for %r:\n%s", makefile["Filename"], error
            )
            return 1

    pipeline = Pypeline(
        nodes=nodes,
        temp_root=config.temp_root,
        max_threads=config.max_threads,
    )

    return pipeline.run(config.pipeline_mode)
