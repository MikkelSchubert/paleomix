#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import time
import logging

import paleomix
import paleomix.logger
import paleomix.resources
import paleomix.yaml

from paleomix.common.console import \
    print_err, \
    print_info

from paleomix.pipeline import \
    Pypeline
from paleomix.nodes.picard import \
    BuildSequenceDictNode
from paleomix.nodes.samtools import \
    FastaIndexNode
from paleomix.nodes.bwa import \
    BWAIndexNode
from paleomix.nodes.bowtie2 import \
    Bowtie2IndexNode
from paleomix.nodes.validation import \
    ValidateFASTAFilesNode

from paleomix.tools.bam_pipeline.makefile import \
    MakefileError, \
    read_makefiles

from paleomix.tools.bam_pipeline.parts import \
    Reads

import paleomix.tools.bam_pipeline.parts as parts
import paleomix.tools.bam_pipeline.config as bam_config
import paleomix.tools.bam_pipeline.mkfile as bam_mkfile


def build_pipeline_trimming(config, makefile):
    """Builds only the nodes required to produce trimmed reads.
    This reduces the required complexity of the makefile to a minimum."""

    nodes = []
    for (_, samples) in makefile["Targets"].iteritems():
        print_info(".", end='')

        for (_, libraries) in samples.iteritems():
            for (_, barcodes) in libraries.iteritems():
                for (barcode, record) in barcodes.iteritems():
                    if record["Type"] in ("Raw", "Trimmed"):
                        offset = record["Options"]["QualityOffset"]
                        reads = Reads(config, record, offset)

                        nodes.extend(reads.nodes)

    return nodes


def build_pipeline_full(config, makefile, return_nodes=True):
    result = []
    features = makefile["Options"]["Features"]
    for (target_name, sample_records) in makefile["Targets"].iteritems():
        print_info(".", end='')

        prefixes = []
        for (_, prefix) in makefile["Prefixes"].iteritems():
            samples = []
            for (sample_name, library_records) in sample_records.iteritems():
                libraries = []
                for (library_name, barcode_records) in library_records.iteritems():
                    lanes = []
                    for (barcode, record) in barcode_records.iteritems():
                        lane = parts.Lane(config, prefix, record, barcode)

                        # ExcludeReads settings may exlude entire lanes
                        if lane.bams:
                            lanes.append(lane)

                    if lanes:
                        libraries.append(parts.Library(config, target_name, prefix, lanes, library_name))

                if libraries:
                    samples.append(parts.Sample(config, prefix, libraries, sample_name))

            if samples:
                prefixes.append(parts.Prefix(config, prefix, samples, features, target_name))

        if prefixes:
            target = parts.Target(config, prefixes, target_name)

            # Construct coverage, depth-histogram, and summary nodes, etc.
            parts.add_statistics_nodes(config, makefile, target)

            if return_nodes:
                # Extra tasks (e.g. coverage, depth-histograms, etc.)
                result.extend(target.nodes)
                # Output BAM files (raw, realigned)
                result.extend(target.bams.itervalues())
            else:
                result.append(target)

    return result


def index_references(config, makefiles):
    references = {}
    references_bwa = {}
    references_bowtie2 = {}
    for makefile in makefiles:
        for subdd in makefile["Prefixes"].itervalues():
            reference = subdd["Reference"]
            if reference not in references:
                # Validation of the FASTA file; not blocking for the other
                # steps, as it is only expected to fail very rarely, but will
                # block subsequent analyses depending on the FASTA.
                valid_node = ValidateFASTAFilesNode(input_files=reference,
                                                    output_file=reference +
                                                    ".validated")
                # Indexing of FASTA file using 'samtools faidx'
                faidx_node = FastaIndexNode(reference)
                # Indexing of FASTA file using 'BuildSequenceDictionary.jar'
                dict_node = BuildSequenceDictNode(config=config,
                                                  reference=reference,
                                                  dependencies=(valid_node,))

                # Indexing of FASTA file using 'bwa index'
                bwa_node = BWAIndexNode(input_file=reference,
                                        dependencies=(valid_node,))
                # Indexing of FASTA file using ''
                bowtie2_node = Bowtie2IndexNode(input_file=reference,
                                                dependencies=(valid_node,))

                references[reference] = (valid_node, faidx_node, dict_node)
                references_bwa[reference] = (valid_node, faidx_node,
                                             dict_node, bwa_node)
                references_bowtie2[reference] = (valid_node, faidx_node,
                                                 dict_node, bowtie2_node)

            subdd["Nodes"] = references[reference]
            subdd["Nodes:BWA"] = references_bwa[reference]
            subdd["Nodes:Bowtie2"] = references_bowtie2[reference]


def run(config, args, pipeline_variant):
    if pipeline_variant not in ("bam", "trim"):
        raise ValueError("Unexpected BAM pipeline variant (%r)"
                         % (pipeline_variant,))

    if not os.path.exists(config.temp_root):
        try:
            os.makedirs(config.temp_root)
        except OSError, error:
            print_err("ERROR: Could not create temp root:\n\t%s" % (error,))
            return 1

    if not os.access(config.temp_root, os.R_OK | os.W_OK | os.X_OK):
        print_err("ERROR: Insufficient permissions for temp root: '%s'"
                  % (config.temp_root,))
        return 1

    # Init worker-threads before reading in any more data
    pipeline = Pypeline(config)

    try:
        print_info("Reading makefiles ...")
        makefiles = read_makefiles(config, args, pipeline_variant)
    except (MakefileError, paleomix.yaml.YAMLError, IOError), error:
        print_err("Error reading makefiles:",
                  "\n  %s:\n   " % (error.__class__.__name__,),
                  "\n    ".join(str(error).split("\n")))
        return 1

    logfile_template = time.strftime("bam_pipeline.%Y%m%d_%H%M%S_%%02i.log")
    paleomix.logger.initialize(config, logfile_template)
    logger = logging.getLogger(__name__)

    pipeline_func = build_pipeline_trimming
    if pipeline_variant == "bam":
        # Build .fai files for reference .fasta files
        index_references(config, makefiles)

        pipeline_func = build_pipeline_full

    print_info("Building BAM pipeline ", end='')
    for makefile in makefiles:
        # If a destination is not specified, save results in same folder as the
        # makefile
        filename = makefile["Statistics"]["Filename"]
        old_destination = config.destination
        if old_destination is None:
            config.destination = os.path.dirname(filename)

        try:
            nodes = pipeline_func(config, makefile)
        except paleomix.node.NodeError, error:
            logger.error("Error while building pipeline for '%s':\n%s",
                         filename, error)
            return 1

        config.destination = old_destination

        pipeline.add_nodes(*nodes)

    print_info("")

    if config.list_input_files:
        logger.info("Printing output files ...")
        pipeline.print_input_files()
        return 0
    elif config.list_output_files:
        logger.info("Printing output files ...")
        pipeline.print_output_files()
        return 0
    elif config.list_executables:
        logger.info("Printing required executables ...")
        pipeline.print_required_executables()
        return 0
    elif config.dot_file:
        logger.info("Writing dependency graph to %r ...", config.dot_file)
        if not pipeline.to_dot(config.dot_file):
            return 1
        return 0

    logger.info("Running BAM pipeline ...")
    if not pipeline.run(dry_run=config.dry_run,
                        max_threads=config.max_threads,
                        progress_ui=config.progress_ui):
        return 1

    return 0


def _print_usage(pipeline):
    basename = "%s_pipeline" % (pipeline,)

    print_info("BAM Pipeline v%s\n" % (paleomix.__version__,))
    print_info("Usage:")
    print_info("  -- %s help           -- Display this message" % basename)
    print_info("  -- %s example [...]  -- Create example project in folder." % basename)
    print_info("  -- %s makefile [...] -- Print makefile template." % basename)
    print_info("  -- %s dryrun [...]   -- Perform dry run of pipeline on provided makefiles." % basename)
    print_info("     %s                   Equivalent to 'bam_pipeline run --dry-run [...]'." % (" " * len(basename),))
    print_info("  -- %s run [...]      -- Run pipeline on provided makefiles." % basename)
    print_info("  -- %s remap [...]    -- Re-map hits from previous alignment." % basename)


def main(argv, pipeline="bam"):
    assert pipeline in ("bam", "trim"), pipeline

    commands = ("makefile", "mkfile", "run",
                "dry_run", "dry-run", "dryrun",
                "remap", "example", "examples")

    if not argv or (argv[0] == "help"):
        _print_usage(pipeline)
        return 0
    elif argv[0] not in commands:
        _print_usage(pipeline)
        return 1
    elif argv[0] in ("mkfile", "makefile"):
        return bam_mkfile.main(argv[1:], pipeline=pipeline)
    elif argv[0] in ("remap", "remap_prefix"):
        # Import here to avoid circular dependency issues
        import paleomix.tools.bam_pipeline.remap as bam_remap

        return bam_remap.main(argv[1:])
    elif argv[0] in ("example", "examples"):
        return paleomix.resources.copy_example("bam_pipeline", argv[1:])

    try:
        config, args = bam_config.parse_config(argv, pipeline)

        if not args[1:]:
            print_err("Please specify at least one makefile!")
            print_err("Use --help for more information.")
            return 1
        elif args and args[0].startswith("dry"):
            config.dry_run = True
    except bam_config.ConfigError, error:
        print_err(error)
        return 1

    return run(config, args[1:], pipeline_variant=pipeline)
