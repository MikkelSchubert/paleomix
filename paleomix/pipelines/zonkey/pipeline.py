#!/usr/bin/python3
#
# Copyright (c) 2016 Mikkel Schubert <MikkelSch@gmail.com>
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
import shutil
import string
import tarfile

import pysam

import paleomix
import paleomix.common.fileutils as fileutils
import paleomix.common.logging
import paleomix.pipelines.bam.mkfile as bam_mkfile
import paleomix.pipelines.zonkey.config as zonkey_config
import paleomix.pipelines.zonkey.database as database
import paleomix.pipelines.zonkey.parts.common as common_nodes
import paleomix.pipelines.zonkey.parts.mitochondria as mitochondria
import paleomix.pipelines.zonkey.parts.nuclear as nuclear
import paleomix.pipelines.zonkey.parts.report as report
import paleomix.pipelines.zonkey.parts.summary as summary
import paleomix.yaml
from paleomix.common.formats.fasta import FASTA
from paleomix.nodes.raxml import RAxMLRapidBSNode
from paleomix.nodes.samtools import BAMIndexNode
from paleomix.pipeline import Pypeline


def run_pipeline(config, nodes, msg):
    pipeline = Pypeline(config)
    pipeline.add_nodes(nodes)

    paleomix.common.logging.initialize(
        log_level=config.log_level, log_file=config.log_file, name="zonkey"
    )

    logger = logging.getLogger(__name__)
    logger.info(msg)

    if config.list_executables:
        pipeline.print_required_executables()
        return True
    elif config.list_output_files:
        pipeline.print_output_files()
        return True
    elif config.list_input_files:
        pipeline.print_input_files()
        return True

    return pipeline.run(max_threads=config.max_threads, dry_run=config.dry_run)


def build_plink_nodes(config, data, root, bamfile, dependencies=()):
    plink = {"root": os.path.join(root, "results", "plink")}

    ped_node = nuclear.BuildTPEDFilesNode(
        output_root=plink["root"],
        table=config.database.filename,
        downsample=config.downsample_to,
        bamfile=bamfile,
        dependencies=dependencies,
    )

    for postfix in ("incl_ts", "excl_ts"):
        plink[postfix] = nuclear.BuildBEDFilesNode(
            output_prefix=os.path.join(plink["root"], postfix),
            tfam=os.path.join(plink["root"], "common.tfam"),
            tped=os.path.join(plink["root"], postfix + ".tped"),
            plink_parameters=config.database.settings["Plink"],
            dependencies=(ped_node,),
        )

    return plink


def build_admixture_nodes(config, data, root, plink):
    nodes = []
    for postfix in ("incl_ts", "excl_ts"):
        bed_node = plink[postfix]

        admix_root = os.path.join(root, "results", "admixture")
        report_root = os.path.join(root, "figures", "admixture")
        for k_groups in sorted(data.groups):
            replicates = []

            input_file = os.path.join(plink["root"], postfix + ".bed")
            for replicate in range(config.admixture_replicates):
                output_root = os.path.join(admix_root, "%02i" % (replicate,))

                node = nuclear.AdmixtureNode(
                    input_file=input_file,
                    output_root=output_root,
                    k_groups=k_groups,
                    groups=data.groups[k_groups],
                    dependencies=(bed_node,),
                )

                replicates.append(node)

            node = nuclear.SelectBestAdmixtureNode(
                replicates=replicates, output_root=admix_root
            )

            if config.admixture_only:
                nodes.append(node)
            else:
                samples = os.path.join(root, "figures", "samples.txt")
                plot = nuclear.AdmixturePlotNode(
                    input_file=os.path.join(
                        admix_root, "%s.%i.Q" % (postfix, k_groups)
                    ),
                    output_prefix=os.path.join(
                        report_root, "%s_k%i" % (postfix, k_groups)
                    ),
                    samples=samples,
                    order=data.sample_order,
                    dependencies=node,
                )

                nodes.append(plot)

    return nodes


def build_treemix_nodes(config, data, root, plink):
    tmix_root = os.path.join(root, "results", "treemix")

    nodes = []
    for postfix in ("incl_ts", "excl_ts"):
        plink_prefix = os.path.join(plink["root"], postfix)
        plink_nodes = plink[postfix]

        freq_node = nuclear.BuildFreqFilesNode(
            output_prefix=plink_prefix,
            input_prefix=os.path.join(plink["root"], postfix),
            tfam=os.path.join(plink["root"], "common.tfam"),
            parameters=config.database.settings["Plink"],
            dependencies=plink_nodes,
        )

        tmix_prefix = os.path.join(tmix_root, postfix)
        tmix_file_node = nuclear.FreqToTreemixNode(
            input_file=plink_prefix + ".frq.strat.gz",
            output_file=tmix_prefix + ".gz",
            dependencies=(freq_node,),
        )

        k_snps = config.treemix_k
        if not k_snps:
            k_snps = (
                "n_sites_%s" % (postfix,),
                os.path.join(plink["root"], "common.summary"),
            )

        for n_migrations in (0, 1):
            n_prefix = "%s.%i" % (tmix_prefix, n_migrations)

            tmix_node = nuclear.TreemixNode(
                data=data,
                input_file=tmix_prefix + ".gz",
                output_prefix=n_prefix,
                m=n_migrations,
                k=k_snps,
                outgroup=config.treemix_outgroup,
                dependencies=(tmix_file_node,),
            )

            samples = os.path.join(root, "figures", "samples.txt")
            output_prefix = os.path.join(
                root, "figures", "treemix", "%s_%i" % (postfix, n_migrations)
            )
            plot_node = nuclear.PlotTreemixNode(
                samples=samples,
                prefix=n_prefix,
                output_prefix=output_prefix,
                dependencies=(tmix_node,),
            )

            nodes.append(plot_node)

    return nodes


def build_pca_nodes(config, data, root, plink):
    pca_root = os.path.join(root, "results", "pca")

    nodes = []
    for postfix in ("incl_ts", "excl_ts"):
        plink_prefix = os.path.join(plink["root"], postfix)
        plink_nodes = plink[postfix]

        pca_prefix = os.path.join(pca_root, postfix)
        pca_node = nuclear.SmartPCANode(
            input_prefix=plink_prefix,
            output_prefix=pca_prefix,
            nchroms=data.settings["NChroms"],
            dependencies=plink_nodes,
        )

        samples = os.path.join(root, "figures", "samples.txt")
        pca_plots = os.path.join(root, "figures", "pca", postfix)
        pca_plot_node = nuclear.PlotPCANode(
            samples=samples,
            prefix=pca_prefix,
            output_prefix=pca_plots,
            dependencies=pca_node,
        )

        nodes.append(pca_plot_node)

    return nodes


def build_coverage_nodes(contigs, mapping, root, nuc_bam, dependencies=()):
    output_prefix = os.path.join(root, "figures", "coverage", "coverage")

    return (
        nuclear.PlotCoverageNode(
            contigs=contigs,
            mapping=mapping,
            input_file=nuc_bam,
            output_prefix=output_prefix,
            dependencies=dependencies,
        ),
    )


def build_mito_nodes(config, root, bamfile, dependencies=()):
    if config.database.mitochondria is None:
        log = logging.getLogger(__name__)
        log.warning("No MT sequences in Zonkey database; cannot perform MT analysis")
        return ()

    samples = os.path.join(root, "figures", "samples.txt")

    mt_prefix = os.path.join(root, "results", "mitochondria", "sequences")
    alignment = mitochondria.MitoConsensusNode(
        database=config.database.filename,
        bamfile=bamfile,
        output_prefix=mt_prefix,
        dependencies=dependencies,
    )

    raxml_template = os.path.join(root, "results", "mitochondria", "raxml_%s")
    phylo = RAxMLRapidBSNode(
        input_alignment=mt_prefix + ".phy",
        output_template=raxml_template,
        model="GTRGAMMA",
        replicates=100,
        dependencies=(alignment,),
    )

    output_prefix = os.path.join(root, "figures", "mitochondria", "mito_phylo")
    trees = mitochondria.DrawPhylogenyNode(
        samples=samples,
        treefile=raxml_template % ("bestTree",),
        bootstraps=raxml_template % ("bootstrap",),
        output_prefix=output_prefix,
        dependencies=(phylo,),
    )

    return (trees,)


def build_pipeline(config, root, nuc_bam, mito_bam, cache):
    nodes = []
    sample_tbl = os.path.join(root, "figures", "samples.txt")
    samples = common_nodes.WriteSampleList(config=config, output_file=sample_tbl)

    if nuc_bam is not None:
        nuc_bam, nuc_bam_info = nuc_bam["Path"], nuc_bam["Info"]

        # When not sampling, BuildTPED relies on indexed access to ease
        # processing of one chromosome at a time. The index is further required
        # for idxstats used by the PlotCoverageNode.
        index = cache.get(nuc_bam)
        if index is None:
            index = cache[nuc_bam] = BAMIndexNode(infile=nuc_bam)

        plink = build_plink_nodes(
            config, config.database, root, nuc_bam, dependencies=(samples, index)
        )

        nodes.extend(build_admixture_nodes(config, config.database, root, plink))

        if not config.admixture_only:
            nodes.extend(
                build_coverage_nodes(
                    contigs=config.database.contigs,
                    mapping=nuc_bam_info.nuclear_contigs,
                    root=root,
                    nuc_bam=nuc_bam,
                    dependencies=(index,),
                )
            )
            nodes.extend(build_pca_nodes(config, config.database, root, plink))
            nodes.extend(build_treemix_nodes(config, config.database, root, plink))

    if mito_bam is not None and not config.admixture_only:
        index = cache.get(mito_bam)
        if index is None:
            index = cache[mito_bam] = BAMIndexNode(infile=mito_bam)

        nodes.extend(
            build_mito_nodes(config, root, mito_bam, dependencies=(samples, index))
        )

    if not config.admixture_only:
        nodes.append(
            report.ReportNode(config, root, nuc_bam, mito_bam, dependencies=nodes)
        )

    return nodes


def run_admix_pipeline(config):
    log = logging.getLogger(__name__)
    log.info("Building %i Zonkey pipeline(s):", len(config.samples))
    config.temp_root = os.path.join(config.destination, "temp")
    if not config.dry_run:
        fileutils.make_dirs(config.temp_root)

    cache = {}
    nodes = []
    items = iter(config.samples.items())
    for idx, (name, sample) in enumerate(sorted(items), start=1):
        root = sample["Root"]
        nuc_bam = sample["Files"].get("Nuc")
        mito_bam = sample["Files"].get("Mito")

        genomes = []
        if mito_bam:
            genomes.append("MT")
        if nuc_bam:
            genomes.append("Nuclear")

        log.info("  %i. %s: %s DNA", idx, name, " and ".join(genomes))

        nodes.extend(build_pipeline(config, root, nuc_bam, mito_bam, cache))

    if config.multisample and not config.admixture_only:
        nodes = [summary.SummaryNode(config, nodes)]

    if not run_pipeline(config, nodes, "Running Zonkey"):
        return 1


def setup_mito_mapping(config):
    genomes_root = os.path.join(config.destination, "genomes")
    if not os.path.exists(genomes_root):
        fileutils.make_dirs(genomes_root)

    mkfile_fpath = os.path.join(config.destination, "makefile.yaml")

    filenames = [mkfile_fpath]
    for name, record in sorted(config.database.mitochondria.items()):
        filenames.append(os.path.join(genomes_root, "%s.fasta" % (record.name,)))

    existing_filenames = [
        filename for filename in filenames if os.path.exists(filename)
    ]

    # A bit strict, but avoid accidential overwrites
    if existing_filenames:
        log = logging.getLogger(__name__)
        log.error("Output file(s) already exists, cannot proceed:")
        for filename in sorted(existing_filenames):
            log.error(" - %r", filename)

        return 1

    with open(mkfile_fpath, "w") as mkfile:
        mkfile.write(
            bam_mkfile.build_makefile(add_prefix_tmpl=False, add_sample_tmpl=False)
        )

        mkfile.write("\n\nPrefixes:\n")

        for name, record in sorted(config.database.mitochondria.items()):
            if "EXCLUDE" in record.meta.upper():
                continue

            mkfile.write("  %s:\n" % (record.name,))
            mkfile.write("    Path: genomes/%s.fasta\n" % (record.name,))

            info = config.database.samples.get(record.name)
            if info is not None:
                mkfile.write("    # Species: %s\n" % (info.get("Species", "NA"),))
                mkfile.write("    # Sex: %s\n" % (info.get("Sex", "NA"),))
                mkfile.write(
                    "    # Publication: %s\n" % (info.get("Publication", "NA"),)
                )
                mkfile.write("    # Sample ID: %s\n" % (info.get("SampleID", "NA"),))

            mkfile.write("\n")

            fasta_fpath = os.path.join(genomes_root, "%s.fasta" % (record.name,))

            with open(fasta_fpath, "w") as fasta_handle:
                record = FASTA(
                    name=record.name,
                    meta=None,
                    sequence=record.sequence.replace("-", ""),
                )

                record.write(fasta_handle)

        mkfile.write("\n")

    return 0


def setup_example(config):
    root = os.path.join(config.destination, "zonkey_pipeline")
    log = logging.getLogger(__name__)
    log.info("Copying example project to %r", root)

    with tarfile.TarFile(config.database.filename) as tar_handle:
        example_files = []
        existing_files = []
        for member in tar_handle.getmembers():
            if os.path.dirname(member.name) == "examples" and member.isfile():
                example_files.append(member)

                destination = fileutils.reroot_path(root, member.name)
                if os.path.exists(destination):
                    existing_files.append(destination)

        if existing_files:
            log.error("Output files already exist at destination:")
            for filename in sorted(existing_files):
                log.error(" - %r", filename)
            return 1
        elif not example_files:
            log.error(
                "Sample database %r does not contain example data; cannot proceed.",
                config.database.filename,
            )
            return 1

        if not os.path.exists(root):
            fileutils.make_dirs(root)

        for member in example_files:
            destination = fileutils.reroot_path(root, member.name)
            src_handle = tar_handle.extractfile(member)
            with open(destination, "wb") as out_handle:
                shutil.copyfileobj(src_handle, out_handle)

    log.info("Sucessfully saved example data in %r", root)

    return 0


def _process_samples(config):
    log = logging.getLogger(__name__)
    for name, info in sorted(config.samples.items()):
        files = {}

        if name == "-":
            log.info("Validating unnamed sample")
        else:
            log.info("Validating sample %r", name)

        for filename in info.pop("Files"):
            filetype = config.database.validate_bam(filename)
            if not filetype:
                log.error("File is not a valid BAM file: %r" % (filename,))
                return False

            if filetype.is_nuclear and filetype.is_mitochondrial:
                if "Nuc" in files:
                    log.error("Two nuclear BAMs specified!")
                    return False
                elif "Mito" in files:
                    log.warning(
                        "Nuclear + mitochondrial BAM and mitochondrial BAM specified; "
                        "the mitochondrial genome in the combined BAM will not be used!"
                    )

                files["Nuc"] = {"Path": filename, "Info": filetype}
                files.setdefault("Mito", filename)
            elif filetype.is_nuclear:
                if "Nuc" in files:
                    log.error("Two nuclear BAMs specified!")
                    return False

                files["Nuc"] = {"Path": filename, "Info": filetype}
            elif filetype.is_mitochondrial:
                if "Mito" in files:
                    log.error("Two nuclear BAMs specified!")
                    return False

                files["Mito"] = filename
            else:
                log.error(
                    "BAM does not contain usable nuclear or mitochondrial contigs: %r",
                    filename,
                )
                return False

        config.samples[name]["Files"] = files

    return True


def _read_sample_table(config, filename):
    """Parses a 2 - 3 column tab-seperated table containing, on each row, a
    name to be used for a sample in the first row, and then the paths two
    either one or to two BAM files, which must represent a single nuclear or
    a single mitochondrial alignment (2 columns), or both (3 columns).
    """
    log = logging.getLogger(__name__)
    log.info("Reading table of samples from %r", filename)
    valid_characters = frozenset(string.ascii_letters + string.digits + ".-_")

    samples = config.samples = {}
    with fileutils.open_ro(filename) as handle:
        for linenum, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            fields = [_f for _f in map(str.strip, line.split("\t")) if _f]
            if len(fields) not in (2, 3):
                log.error(
                    "Error reading sample table (%r) at line %i: Expected 2 or 3 "
                    "columns, found %i; please correct file before continuing.",
                    filename,
                    linenum,
                    len(fields),
                )
                return

            name = fields[0]
            invalid_letters = frozenset(name) - valid_characters
            if invalid_letters:
                log.error(
                    "Error reading sample table (%r) at line %i: Sample name contains "
                    "illegal character(s). Only letters, numbers, and '-', '_', and "
                    "'.' are allowed, but found %r in name %r ",
                    filename,
                    linenum,
                    "".join(invalid_letters),
                    name,
                )
                return
            elif name in samples:
                log.error(
                    "Duplicate name %r in sample table; names must be unique!", name
                )
                return

            samples[name] = {
                "Root": os.path.join(config.destination, name),
                "Files": fields[1:],
            }

    return True


def finalize_run_config(parser, args):
    log = logging.getLogger(__name__)
    if args.command in ("run", "dryrun") and not (1 <= len(args.files) <= 3):
        parser.print_usage()
        return

    args.dry_run |= args.command == "dryrun"
    args.multisample = False

    known_samples = set(args.database.samples) | set(("Sample",))
    unknown_samples = set(args.treemix_outgroup) - known_samples
    if unknown_samples:
        log.error(
            "Argument --treemix-outgroup includes unknown sample(s): %s; known "
            "samples are %s. Note that names are case-sensitive."
            ", ".join(map(repr, sorted(unknown_samples))),
            ", ".join(map(repr, sorted(known_samples))),
        )
        return

    if len(args.files) == 1:
        args.files.append(fileutils.swap_ext(args.files[0], ".zonkey"))

    if len(args.files) == 2:
        filename, args.destination = args.files

        if os.path.exists(args.destination) and not os.path.isdir(args.destination):
            log.error("Destination %r is not a directory", args.destination)
            return
        elif not os.path.isfile(filename):
            log.error("Not a valid filename: %r", filename)
            return
        elif _is_bamfile(filename):
            args.samples = {"-": {"Root": args.destination, "Files": [filename]}}
        else:
            args.multisample = True
            if not _read_sample_table(args, filename):
                return
    elif len(args.files) == 3:
        filename_1, filename_2, args.destination = args.files

        args.samples = {
            "-": {"Root": args.destination, "Files": [filename_1, filename_2]}
        }
    else:
        raise RuntimeError("Unexpected number of arguments: %r" % (args.files,))

    # Identify (mito or nuc?) and validate BAM files provided by user
    if not _process_samples(args):
        return

    return args


def _is_bamfile(filename):
    """Returns true if a file is a BAM file, false otherwise.
    """
    try:
        with pysam.AlignmentFile(filename, "rb"):
            return True
    except ValueError:
        return False
    except IOError:
        return False


def main(argv):
    parser, run_parser = zonkey_config.build_parser()
    if not argv:
        parser.print_help()
        return

    args = parser.parse_args(argv)
    log = logging.getLogger(__name__)

    try:
        args.database = database.ZonkeyDB(args.database)
    except database.ZonkeyDBError as error:
        log.error("Error reading database %r: %s", args.database, error)
        return

    if args.command in ("run", "dryrun"):
        args = finalize_run_config(run_parser, args)
        if args is not None:
            return run_admix_pipeline(args)
    elif args.command == "mito":
        return setup_mito_mapping(args)
    elif args.command == "example":
        return setup_example(args)

    return 1
