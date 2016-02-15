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
"""
This is a tool to ease the process of submitting data previously processed.
using the BAM Pipeline to ENA. This accomplished in roughly 5 steps:

1. Create the project on ENA and create each of the samples that are to be
   uploaded for that project.

2. Run 'paleomix ena build-table' on one or more makefiles; this generates a
   basic template for the lanes to be uploaded. This table is edited to
   reflect the desired layout of uploaded data (except the MAKEFILE_PATH
   column, which MUST NOT be edited).

3. Run 'paleomix ena collect-files' on the table built in step 2; this collects
   all FASTQ files into a single folder, merging lanes into one (for SE) or two
   (for PE) files and calculates MD5 sums. The folder containing the merged
   FASTQ files is uploaded to ENA, using the same path (e.g. upload to
   /ena_submission/fastq/ by default).

4. Run 'paleomix ena finalize-table' on the table, with the 'se' and 'pe'
   options. This is required since the two types of tables differ in layout.

5. Uploaded the resulting 'se' and 'pe' tables to ENA.
"""
import argparse
import glob
import os
import string
import sys

import paleomix.yaml

import paleomix.ui as ui

import paleomix.tools.factory as factory

from paleomix.atomiccmd.command import \
    AtomicCmd

from paleomix.atomiccmd.sets import \
    ParallelCmds

from paleomix.node import \
    CommandNode

from paleomix.pipeline import \
    Pypeline

import paleomix.common.fileutils as fileutils


###############################################################################

class ENAError(RuntimeError):
    pass


###############################################################################

class CatFilesNode(CommandNode):
    """Runs the equivalent of 'paleomix cat $@ | gzip > ${DST}'."""
    def __init__(self, input_files, destination, dependencies=()):
        cat_cmd = factory.new("cat")
        cat_cmd.add_multiple_values(input_files)
        cat_cmd.set_kwargs(OUT_STDOUT=AtomicCmd.PIPE)
        cat_cmd = cat_cmd.finalize()

        zip_cmd = AtomicCmd("gzip",
                            IN_STDIN=cat_cmd,
                            OUT_STDOUT=destination)

        description = "<Cat %s -> %s>" \
            % (fileutils.describe_files(input_files), destination)

        CommandNode.__init__(self,
                             description=description,
                             command=ParallelCmds((cat_cmd, zip_cmd)),
                             dependencies=dependencies)


class MD5FilesNode(CommandNode):
    """Runs the equivalent of 'md5sum $1 > $2'."""
    def __init__(self, input_file, destination, dependencies=()):
        md5_cmd = AtomicCmd(("md5sum", "%(IN_FILE)s"),
                            IN_FILE=input_file,
                            OUT_STDOUT=destination)

        description = "<MD5Sum %s -> %s>" \
            % (input_file, destination)

        CommandNode.__init__(self,
                             description=description,
                             command=md5_cmd,
                             dependencies=dependencies)


###############################################################################

def _is_paired_end(template):
    """Returns true if a path from a makefile is for paired-end data."""
    return (template.format(Pair=1) != template)


def _sorted_glob(tmpl):
    return list(sorted(glob.glob(tmpl))) or [tmpl]


def _collect_files(template):
    if not _is_paired_end(template):
        return {"": _sorted_glob(template)}

    files = {}
    for (mate, name) in enumerate(("_R1", "_R2"), start=1):
        files[name] = _sorted_glob(template.format(Pair=mate))

    assert len(files["_R1"]) == len(files["_R2"]), template

    return files


def _build_filename(args, row, extensions, postfix=""):
    return os.path.join(args.destination,
                        extensions.split(".")[0],
                        row["MAKEFILE_TARGET"],
                        "%s__%s__%s%s.%s"
                        % (row["MAKEFILE_SAMPLE"],
                           row["MAKEFILE_LIBRARY"],
                           row["MAKEFILE_LANE"],
                           postfix,
                           extensions))


###############################################################################


def _main_collect_files(args):
    """Main function for collecting FASTQ files and calculating MD5 sums."""
    args.temp_root = os.path.join(args.destination, "temp")

    nodes = _build_collect_nodes(args)
    if nodes is None:
        return 1

    pipeline = Pypeline(args)
    pipeline.add_nodes(nodes)

    if args.list_output_files:
        pipeline.print_output_files()
    elif args.list_executables:
        pipeline.print_required_executables()
    elif not pipeline.run(max_threads=args.max_threads,
                          dry_run=args.dry_run):
        return 1

    return 0


def _build_collect_nodes(args):
    """Builds a set of nodes for merging and gzip-compressing the reads from
    each lane, as well as calculating the MD5 sums for each resulting file. By
    default, the resulting files are located in the folders
        - ena_submission/fastq
        - ena_submission/md5
    """
    nodes = []

    for row in _read_table(args.table):
        path = row["MAKEFILE_PATH"]

        for postfix, filenames in sorted(_collect_files(path).iteritems()):
            destination_fastq = _build_filename(args, row, "fastq.gz", postfix)
            destination_md5 = _build_filename(args, row, "md5", postfix)

            for filename in filenames:
                if not os.path.exists(filename):
                    _build_collect_nodes_error(row, "File not found", filename)
                elif not os.path.isfile(filename):
                    _build_collect_nodes_error(row, "Not a file", filename)

            cat_node = CatFilesNode(input_files=filenames,
                                    destination=destination_fastq)

            md5_node = MD5FilesNode(input_file=destination_fastq,
                                    destination=destination_md5,
                                    dependencies=(cat_node,))

            nodes.append(md5_node)

    return nodes


def _build_collect_nodes_error(row, msg, filename):
    row = dict(row)
    row["PATH"] = filename
    row["MSG"] = msg

    raise ENAError(("%(MSG)s; "
                    "target = %(MAKEFILE_TARGET)s; "
                    "sample = %(MAKEFILE_SAMPLE)s; "
                    "library = %(MAKEFILE_LIBRARY)s; "
                    "lane = %(MAKEFILE_LANE)s; "
                    "path = %(PATH)s") % row)


###############################################################################

def _main_build_table(args):
    keys = ("MAKEFILE_TARGET", "MAKEFILE_SAMPLE", "MAKEFILE_LIBRARY",
            "MAKEFILE_LANE", "sample_alias", "instrument_model",
            "library_source", "library_selection",
            "library_strategy", "design_description",
            "library_construction_protocol", "insert_size",
            "MAKEFILE_PATH")

    result = []
    for row in _build_table_rows(args):
        result.append([row[key] for key in keys])

    print "\t".join(keys)
    for row in sorted(result):
        print "\t".join(row)


def _build_table_rows(args):
    rows = []
    for filename in args.makefile:
        for (target, sample, library, lane, path) in _parse_makefile(filename):
            if isinstance(path, dict):
                ui.print_err("WARNING: Found pre-processed data "
                             "at %s:%s:%s:%s; cannot collect raw "
                             "FASTQ data."
                             % (target, sample, library, lane))
                continue

            row = {"sample_alias": "*",
                   "instrument_model": "*",
                   "library_source": "GENOMIC",
                   "library_selection": "RANDOM",
                   "library_strategy": "WGS",
                   "design_description": "",
                   "library_construction_protocol": "",
                   "insert_size": "0",

                   "MAKEFILE_TARGET": target,
                   "MAKEFILE_SAMPLE": sample,
                   "MAKEFILE_LIBRARY": library,
                   "MAKEFILE_LANE": lane,
                   "MAKEFILE_PATH": path}

            rows.append(row)

    return rows


def _parse_makefile(filename):
    with open(filename) as handle:
        mkfile = paleomix.yaml.safe_load(handle.read())
        mkfile.pop("Options", None)
        mkfile.pop("Prefixes", None)

        for target, samples in sorted(mkfile.iteritems()):
            samples.pop("Options", None)

            for sample, libraries in sorted(samples.iteritems()):
                libraries.pop("Options", None)

                for library, lanes in sorted(libraries.iteritems()):
                    lanes.pop("Options", None)

                    for lane, path in sorted(lanes.iteritems()):
                        yield (target, sample, library, lane, path)


###############################################################################

def _main_finalize_table(args):
    table = _read_table(args.table)
    for row in table:
        if _is_paired_end(row["MAKEFILE_PATH"]):
            if args.mode == "pe":
                row['forward_file_name'] = _build_filename(args, row,
                                                           "fastq.gz", "_R1")
                row['forward_file_md5'] = read_md5(args, row, "_R1")
                row['reverse_file_name'] = _build_filename(args, row,
                                                           "fastq.gz", "_R2")
                row['reverse_file_md5'] = read_md5(args, row, "_R2")
        elif args.mode == "se":
            row["file_name"] = _build_filename(args, row, "fastq.gz")
            row["file_md5"] = read_md5(args, row)

        row.setdefault("library_name", row["MAKEFILE_LIBRARY"])

    keys = ("sample_alias", "instrument_model", "library_name",
            "library_source", "library_selection", "library_strategy",
            "design_description", "library_construction_protocol")

    if args.mode == "pe":
        keys += ('insert_size',
                 'forward_file_name', 'forward_file_md5',
                 'reverse_file_name', 'reverse_file_md5')
    else:
        keys += ('file_name', 'file_md5')

    results = []
    do_warn = False
    for row in table:
        if _is_paired_end(row["MAKEFILE_PATH"]) == (args.mode == "pe"):
            results.append([row[key] for key in keys])
        else:
            do_warn = True

    print "\t".join(keys)
    for row in sorted(results):
        print "\t".join(row)

    if do_warn:
        current, other = "single-end", "paired-end"
        if args.mode == "pe":
            other, current = current, other

        ui.print_warn("WARNING: Table for {current} has been printed, but\n"
                      "         the table also contains data for {other}\n"
                      "         reads. Remember to construct a table for\n"
                      "         {other} reads as well!"
                      .format(current=current, other=other))

    return 0


def read_md5(args, row, postfix=""):
    """Reads the md5-sum generated for a given lane; arguments correspond to
    'build_filename', excepting that the extension is already defined.
    """
    filename = _build_filename(args, row, "md5", postfix)
    with open(filename) as handle:
        lines = handle.readlines()

        if len(lines) != 1:
            raise ENAError("MD5-sum file (%r) does not match expected format; "
                           "exactly one line expected, but found %i line(s); "
                           "please remove this file and re-run the "
                           "'collect-files' command!"
                           % (filename, len(lines)))

        fields = lines[0].split()
        if len(fields) != 2:
            raise ENAError("MD5-sum file (%r) does not match expected format; "
                           "exactly two columns expected, but found %i; "
                           "please remove this file and re-run the "
                           "'collect-files' command!"
                           % (filename, len(fields)))

        if len(fields[0]) != 32 or (set(fields[0]) - set(string.hexdigits)):
            raise ENAError("MD5-sum file (%r) does not match expected format; "
                           "32 digit hexadecimal expected, but found %r; "
                           "please remove this file and re-run the "
                           "'collect-files' command!"
                           % (filename, fields[0]))

        return fields[0]


def _read_table(filename):
    """Reads table generated using 'build-table'."""
    rows = []
    with open(filename) as handle:
        header = None
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if header is None:
                header = fields
            elif fields and line.strip():
                if len(fields) == len(header):
                    rows.append(dict(zip(header, fields)))
                else:
                    assert False

    return rows


###############################################################################

def parse_args(argv):
    """Parses the result of sys.argv[1:] or equivalent.

    The resulting object contains a 'function' value, corresponding to the
    command supplied by the user. This function should simply be called with
    the parse_args result as its only parameter.
    """

    description = "This is a tool to ease the process of submitting data " \
        "previously processed using the BAM Pipeline to ENA. Please see the " \
        "documentation for instructions on how to use."

    parser = argparse.ArgumentParser(prog="paleomix ena",
                                     description=description)

    subparsers = parser.add_subparsers()

    ###########################################################################
    # Parser for the 'build-table' command
    parser_build_table = subparsers.add_parser('build-table')
    parser_build_table.add_argument('makefile', nargs="+",
                                    help="One or more BAM Pipeline makefiles.")
    parser_build_table.set_defaults(function=_main_build_table)

    ###########################################################################
    # Parser for the 'collect-files' command
    parser_collect = subparsers.add_parser('collect-files')
    parser_collect.add_argument('table',
                                help="Table generated using the 'build-table' "
                                     "command.")
    parser_collect.add_argument('--destination', default="ena_submission",
                                metavar="FOLDER",
                                help="Destination folder for temporary files, "
                                     "merged FASTQ files, and MD5-sum files "
                                     "[default: %(default)s].")
    parser_collect.add_argument('--max-threads', default=2, type=int,
                                metavar="N",
                                help="Maximum number of simultanous steps "
                                     "to execute [default: %(default)s].")
    parser_collect.add_argument('--list-output-files', default=False,
                                action="store_true",
                                help="List the (status of) files (to be) "
                                     "generated by this command.")
    parser_collect.add_argument('--list-executables', default=False,
                                action="store_true",
                                help="List all executables required by the "
                                     "pipeline, with version requirements "
                                     "(if any).")
    parser_collect.add_argument("--dry-run",
                                action="store_true", default=False,
                                help="If passed, only a dry-run in performed, "
                                     "the dependency tree is printed, and no "
                                     "tasks are executed.")
    parser_collect.set_defaults(function=_main_collect_files)

    ###########################################################################
    # Parser for the 'finalize-table' command
    parser_finalize_table = subparsers.add_parser('finalize-table')
    parser_finalize_table.add_argument('mode', choices=("se", "pe"),
                                       help="Output table containing either "
                                            "the single-end ('se') or paired"
                                            "-end ('pe') reads in the input "
                                            "table.")
    parser_finalize_table.add_argument('table',
                                       help="Table generated using the "
                                            "'build-table' command.")
    parser_finalize_table.add_argument('--destination',
                                       default="ena_submission",
                                       metavar="FOLDER",
                                       help="Destination folder for temporary "
                                            "files, merged FASTQ files, and "
                                            "MD5-sum files "
                                            "[default: %(default)s/].")
    parser_finalize_table.set_defaults(function=_main_finalize_table)

    return parser.parse_args(argv)


def main(argv):
    """Main function; takes a list of arguments but excluding sys.argv[0]."""
    args = parse_args(argv)

    try:
        return args.function(args)
    except ENAError, error:
        ui.print_err("FATAL ERROR:\n  %s" % (error,))

    return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
