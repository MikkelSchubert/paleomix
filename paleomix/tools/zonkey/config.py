#!/usr/bin/python
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
import optparse
import os
import string
import sys

import pysam

import paleomix
import paleomix.common.fileutils as fileutils
import paleomix.tools.zonkey.database as database
from paleomix.config import PerHostConfig, PerHostValue

_USAGE = """USAGE:
{0} run <SampleDB> <samples.txt> [<destination>]
{0} run <SampleDB> <sample.bam> [<destination>]
{0} run <SampleDB> <nuclear.bam> <mitochondrial.bam> <destination>
{0} dryrun <SampleDB> [...]
{0} mito <SampleDB> <destination>
{0} example <SampleDB> <destination>
"""

# List of valid commands / aliases, currently mostly undocumented
_CMD_ALIASES = {
    "mito": "mito",
    "mitochondria": "mito",
    "mitochondrial": "mito",
    "nuc": "run",
    "nuclear": "run",
    "run": "run",
    "dryrun": "dryrun",
    "dry_run": "dryrun",
    "example": "example",
    "examples": "example",
}


def print_usage(out=sys.stderr):
    out.write(_USAGE.format("paleomix zonkey"))


def parse_config(argv):
    config, args = _parse_arguments(argv)
    if not args:
        print_usage()
        return

    config.command = _CMD_ALIASES.get(args[0])
    if config.command is None:
        log = logging.getLogger(__name__)
        log.error("Unknown command %r" % (args[0],))
        return
    elif config.command == "dryrun":
        config.command = "run"
        config.dry_run = True

    return parse_run_config(config, args[1:])


def parse_run_config(config, args):
    if not (2 <= len(args) <= 4):
        print_usage()
        return

    log = logging.getLogger(__name__)
    config.multisample = False
    config.tablefile = args[0]

    try:
        config.database = database.ZonkeyDB(config.tablefile)
    except database.ZonkeyDBError as error:
        log.error("Error reading database %r: %s", config.tablefile, error)
        return

    known_samples = set(config.database.samples) | set(("Sample",))
    unknown_samples = set(config.treemix_outgroup) - known_samples
    if unknown_samples:
        log.error(
            "Argument --treemix-outgroup includes unknown sample(s): %s; known "
            "samples are %s. Note that names are case-sensitive."
            ", ".join(map(repr, sorted(unknown_samples))),
            ", ".join(map(repr, sorted(known_samples))),
        )
        return

    if config.command in ("mito", "example"):
        if len(args) != 2:
            log.error("Wrong number of arguments")
            print_usage()
            return

        config.destination = args[1]
        config.samples = {}
    elif len(args) == 2:
        filename = args[1]
        config.destination = fileutils.swap_ext(filename, ".zonkey")

        if not os.path.isfile(filename):
            log.error("Not a valid filename: %r", filename)
            return
        elif _is_bamfile(filename):
            # Called as either of
            #   zonkey run <SampleDB> <nuclear.bam>
            #   zonkey run <SampleDB> <mitochondrial.bam>
            config.samples = {"-": {"Root": config.destination, "Files": [filename]}}
        else:
            config.multisample = True
            if not _read_sample_table(config, filename):
                return
    elif 3 <= len(args) <= 4:
        root = args[-1]
        if os.path.exists(root) and not os.path.isdir(root):
            log.error("Missing destination folder.")
            print_usage()
            return

        config.destination = root

        if len(args) == 3:
            # zonkey run <SampleDB> <nuclear.bam> <destination>
            # zonkey run <SampleDB> <mitochondrial.bam> <destination>
            # zonkey run <SampleDB> <samples.txt> <destination>
            filename = args[-2]

            if not os.path.isfile(filename):
                log.error("Not a valid filename: %r", filename)
                return
            elif _is_bamfile(filename):
                # Called as either of
                #   zonkey run <SampleDB> <nuclear.bam>
                #   zonkey run <SampleDB> <mitochondrial.bam>
                config.samples = {
                    "-": {"Root": config.destination, "Files": [filename]}
                }
            else:
                config.multisample = True
                if not _read_sample_table(config, filename):
                    return
        else:
            # zonkey run <SampleDB> <nuclear.bam> <mitochondrial.bam> <dst>
            config.destination = root
            config.samples = {"-": {"Root": root, "Files": args[1:-1]}}
    else:
        raise RuntimeError(
            "Unhandled number of args in parse_config: %i\n" % (len(args),)
        )

    # Identify (mito or nuc?) and validate BAM files provided by user
    if not _process_samples(config):
        return

    return config


def _is_bamfile(filename):
    """Returns true if a file is a BAM file, false otherwise.
    """
    try:
        with pysam.Samfile(filename, "rb"):
            return True
    except ValueError:
        return False
    except IOError:
        return False


def _process_samples(config):
    log = logging.getLogger(__name__)
    for name, info in sorted(config.samples.items()):
        files = {}

        if name == "-":
            log.info("Validating unnamed sample ...")
        else:
            log.info("Validating sample %r ...", name)

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

                files["Nuc"] = filename
                files.setdefault("Mito", filename)
            elif filetype.is_nuclear:
                if "Nuc" in files:
                    log.error("Two nuclear BAMs specified!")
                    return False

                files["Nuc"] = filename
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
    valid_characters = frozenset(string.letters + string.digits + ".-_")

    samples = config.samples = {}
    with fileutils.open_ro(filename) as handle:
        for linenum, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            fields = [_f for _f in map(str.strip, line.split("\t")) if _f]
            if len(fields) not in (2, 3):
                log.error(
                    "Error reading sample table (%r) at line %i:  Expected 2 or 3 "
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


def _parse_arguments(argv):
    per_host_cfg = PerHostConfig("zonkey")

    usage_str = "%prog <command> [options] <SampleDB> <bam/sam> [destination]"
    version_str = "%%prog v%s" % (paleomix.__version__,)
    parser = optparse.OptionParser(usage=usage_str, version=version_str)

    group = optparse.OptionGroup(parser, "Program options")
    group.add_option(
        "--downsample-to",
        type=int,
        default=PerHostValue(1000000),
        help="Number of reads to use for analyses; if 0, no "
        "downsampling is performed [default: %default]",
    )
    group.add_option(
        "--admixture-replicates",
        type=int,
        default=PerHostValue(1),
        help="Number of admixture replicates to run, before "
        "the result with the highest likelihood [%default]",
    )
    group.add_option(
        "--treemix-k",
        type=int,
        default=PerHostValue(0),
        help="Value passed to treemix's -k option; number of "
        "SNPs per block for estimation of the covariance "
        "matrix. If set to 0, a value will be estimated "
        "assuming an even distribution of SNPs [%default]",
    )
    group.add_option(
        "--treemix-outgroup",
        default="",
        help="Comma-seperated list of samples to use as the "
        "outgroup when running TreeMix; note that these "
        "must form a monophyletic clade, or TreeMix will "
        "not execute.",
    )

    group.add_option(
        "--admixture-only",
        help=optparse.SUPPRESS_HELP,
        default=False,
        action="store_true",
    )
    parser.add_option_group(group)

    paleomix.logger.add_optiongroup(parser, default=PerHostValue("warning"))

    group = optparse.OptionGroup(parser, "Pipeline")
    group.add_option(
        "--dry-run",
        action="store_true",
        default=False,
        help="If passed, only a dry-run in performed, the "
        "dependency tree is printed, and no tasks are "
        "executed.",
    )
    group.add_option(
        "--max-threads",
        type=int,
        default=PerHostValue(1),
        help="Maximum number of threads to use [%default]",
    )
    group.add_option(
        "--list-input-files",
        action="store_true",
        default=False,
        help="List all input files used by pipeline for the "
        "makefile(s), excluding any generated by the "
        "pipeline itself.",
    )
    group.add_option(
        "--list-output-files",
        action="store_true",
        default=False,
        help="List all output files generated by pipeline for " "the makefile(s).",
    )
    group.add_option(
        "--list-executables",
        action="store_true",
        default=False,
        help="List all executables required by the pipeline, "
        "with version requirements (if any).",
    )
    parser.add_option_group(group)

    config, args = per_host_cfg.parse_args(parser, argv)

    config.treemix_outgroup = tuple(
        [_f for _f in sorted(config.treemix_outgroup.split(",")) if _f]
    )

    return config, args
