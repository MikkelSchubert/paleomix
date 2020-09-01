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


import os
import logging

import paleomix
import paleomix.common.logging
import paleomix.resources
import paleomix.yaml

from paleomix.pipeline import Pypeline
from paleomix.nodes.samtools import FastaIndexNode
from paleomix.nodes.bwa import BWAIndexNode
from paleomix.nodes.bowtie2 import Bowtie2IndexNode
from paleomix.nodes.validation import ValidateFASTAFilesNode

from paleomix.pipelines.bam.makefile import MakefileError, read_makefiles

from paleomix.pipelines.bam.parts import Reads

import paleomix.pipelines.bam.parts as parts
import paleomix.pipelines.bam.config as bam_config
import paleomix.pipelines.bam.mkfile as bam_mkfile


def build_pipeline_trimming(config, makefile):
    """Builds only the nodes required to produce trimmed reads.
    This reduces the required complexity of the makefile to a minimum."""

    nodes = []
    for (_, samples) in makefile["Targets"].items():
        for libraries in samples.values():
            for barcodes in libraries.values():
                for record in barcodes.values():
                    if record["Type"] in ("Raw", "Trimmed"):
                        offset = record["Options"]["QualityOffset"]
                        reads = Reads(config, record, offset)

                        nodes.extend(reads.nodes)

    return nodes


def build_pipeline_full(config, makefile, return_nodes=True):
    result = []
    features = makefile["Options"]["Features"]
    for (target_name, sample_records) in makefile["Targets"].items():
        prefixes = []
        for (_, prefix) in makefile["Prefixes"].items():
            samples = []
            for (sample_name, library_records) in sample_records.items():
                libraries = []
                for (library_name, barcode_records) in library_records.items():
                    lanes = []
                    for (barcode, record) in barcode_records.items():
                        lane = parts.Lane(config, prefix, record, barcode)

                        # ExcludeReads settings may exlude entire lanes
                        if lane.bams:
                            lanes.append(lane)

                    if lanes:
                        libraries.append(
                            parts.Library(
                                config=config,
                                target=target_name,
                                prefix=prefix,
                                lanes=lanes,
                                name=library_name,
                            )
                        )

                if libraries:
                    samples.append(
                        parts.Sample(
                            config=config,
                            prefix=prefix,
                            libraries=libraries,
                            name=sample_name,
                        )
                    )

            if samples:
                prefixes.append(
                    parts.Prefix(
                        config=config,
                        prefix=prefix,
                        samples=samples,
                        features=features,
                        target=target_name,
                    )
                )

        if prefixes:
            target = parts.Target(config, prefixes, target_name)

            # Construct coverage, depth-histogram, and summary nodes, etc.
            parts.add_statistics_nodes(config, makefile, target)

            if return_nodes:
                # Extra tasks (e.g. coverage, depth-histograms, etc.)
                result.extend(target.nodes)
                # Output BAM files
                result.extend(target.bams.values())
            else:
                result.append(target)

    return result


def index_references(config, makefiles):
    references = {}
    references_bwa = {}
    references_bowtie2 = {}
    for makefile in makefiles:
        for subdd in makefile["Prefixes"].values():
            reference = subdd["Reference"]
            if reference not in references:
                # Validation of the FASTA file; not blocking for the other
                # steps, as it is only expected to fail very rarely, but will
                # block subsequent analyses depending on the FASTA.
                valid_node = ValidateFASTAFilesNode(
                    input_files=reference, output_file=reference + ".validated"
                )
                # Indexing of FASTA file using 'samtools faidx'
                faidx_node = FastaIndexNode(reference)

                # Indexing of FASTA file using 'bwa index'
                bwa_node = BWAIndexNode(
                    input_file=reference, dependencies=(valid_node,)
                )
                # Indexing of FASTA file using ''
                bowtie2_node = Bowtie2IndexNode(
                    input_file=reference, dependencies=(valid_node,)
                )

                references[reference] = (valid_node, faidx_node)
                references_bwa[reference] = (valid_node, faidx_node, bwa_node)
                references_bowtie2[reference] = (valid_node, faidx_node, bowtie2_node)

            subdd["Nodes"] = references[reference]
            subdd["Nodes:BWA"] = references_bwa[reference]
            subdd["Nodes:Bowtie2"] = references_bowtie2[reference]


def run(config, pipeline_variant):
    paleomix.common.logging.initialize(
        log_level=config.log_level, log_file=config.log_file, name="bam_pipeline"
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

    # Init worker-threads before reading in any more data
    pipeline = Pypeline(config)

    try:
        makefiles = read_makefiles(config.makefiles, pipeline_variant)
    except (MakefileError, paleomix.yaml.YAMLError, IOError) as error:
        logger.error("Error reading makefiles: %s", error)
        return 1

    pipeline_func = build_pipeline_trimming
    if pipeline_variant == "bam":
        # Build .fai files for reference .fasta files
        index_references(config, makefiles)

        pipeline_func = build_pipeline_full

    for makefile in makefiles:
        logger.info("Building BAM pipeline for %r", makefile["Filename"])
        try:
            nodes = pipeline_func(config, makefile)
        except paleomix.node.NodeError as error:
            logger.error(
                "Error while building pipeline for %r:\n%s", makefile["Filename"], error
            )
            return 1

        pipeline.add_nodes(*nodes)

    if config.list_input_files:
        logger.info("Printing output files")
        pipeline.print_input_files()
        return 0
    elif config.list_output_files:
        logger.info("Printing output files")
        pipeline.print_output_files()
        return 0
    elif config.list_executables:
        logger.info("Printing required executables")
        pipeline.print_required_executables()
        return 0

    logger.info("Running BAM pipeline")
    if not pipeline.run(dry_run=config.dry_run, max_threads=config.max_threads):
        return 1

    return 0


def main(argv, pipeline="bam"):
    if pipeline not in ("bam", "trim"):
        raise ValueError(pipeline)

    parser = bam_config.build_parser(pipeline)
    if not argv:
        parser.print_help()
        return 0

    args = parser.parse_args(argv)
    if args.command in ("makefile", "mkfile"):
        return bam_mkfile.main(args, pipeline=pipeline)
    elif args.command in ("example",):
        return paleomix.resources.copy_example("bam_pipeline", args)

    if args.command.startswith("dry"):
        args.dry_run = True

    return run(args, pipeline_variant=pipeline)
