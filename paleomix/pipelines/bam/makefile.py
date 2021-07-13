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
import collections
import copy
import glob
import itertools
import logging
import os
import string

import paleomix.pipelines.bam.paths as paths

from paleomix.common.bamfiles import BAM_PLATFORMS
from paleomix.common.fileutils import get_files_glob
from paleomix.common.utilities import fill_dict
from paleomix.common.makefile import (
    MakefileError,
    REQUIRED_VALUE,
    WithoutDefaults,
    read_makefile,
    IsInt,
    IsUnsignedInt,
    IsFloat,
    IsStr,
    IsNone,
    IsBoolean,
    And,
    Or,
    Not,
    ValueIn,
    ValuesIntersect,
    ValuesSubsetOf,
    StringIn,
    StringStartsWith,
    IsListOf,
    IsDictOf,
)
from paleomix.common.formats.fasta import FASTA, FASTAError

import paleomix.common.sequences as sequences


_READ_TYPES = set(("Single", "Singleton", "Collapsed", "CollapsedTruncated", "Paired"))

# The maximum reference sequence length supported by the BAI index format:
#   https://samtools.github.io/hts-specs/SAMv1.pdf
_BAM_MAX_SEQUENCE_LENGTH = 2 ** 29 - 1


def read_makefiles(filenames, pipeline_variant="bam"):
    if pipeline_variant not in ("bam", "trim"):
        raise ValueError(
            "'pipeline_variant' must be 'bam' or 'trim', not %r" % (pipeline_variant,)
        )

    logger = logging.getLogger(__name__)

    makefiles = []
    for filename in filenames:
        logger.info("Reading makefile %r", filename)
        makefile = read_makefile(filename, _VALIDATION)
        makefile = _mangle_makefile(makefile, pipeline_variant)
        makefile["Filename"] = filename

        makefiles.append(makefile)

    return _validate_makefiles(makefiles)


def _alphanum_check(whitelist, min_len=1):
    description = "characters a-z, A-Z, 0-9%s allowed"
    description %= (", and %r" % whitelist,) if whitelist else ""

    whitelist += string.ascii_letters + string.digits

    return And(
        IsStr(min_len=min_len), ValuesSubsetOf(whitelist, description=description)
    )


# Valid names for genomes
_VALID_GENOME_NAME = And(
    _alphanum_check(whitelist="._-*"),
    Not(StringIn(["Options"] + [(s + "Reads") for s in _READ_TYPES])),
)

# Valid paths for genomes; avoids some problems with e.g. Bowtie2
_VALID_GENOME_PATH = And(
    IsStr(), Not(ValuesIntersect('\\:?"<>|() \t\n\v\f\r')), default=REQUIRED_VALUE
)

# Valid strings for samples / libraries / lanes
_VALID_FILENAME = _alphanum_check(whitelist="._-", min_len=2)

# Features that can be specified at project, sample, and library level
_VALID_LIBRARY_LEVEL_FEATURES = {
    "mapDamage",
    "PCRDuplicates",
}

_VALID_FEATURES_DICT = {
    "Coverage": IsBoolean(default=True),
    "Depths": IsBoolean(default=True),
    "mapDamage": StringIn(("rescale", "model", "plot", True, False), default=False),
    "PCRDuplicates": StringIn((True, False, "mark", "filter"), default="filter"),
    "Summary": IsBoolean(default=True),
}

_VALID_EXCLUDE_DICT = {
    "Single": IsBoolean(default=False),
    "Collapsed": IsBoolean(default=False),
    "CollapsedTruncated": IsBoolean(default=False),
    "Paired": IsBoolean(default=False),
    "Singleton": IsBoolean(default=False),
}


_VALIDATION_OPTIONS = {
    # Sequencing platform, used to tag read-groups.
    "Platform": StringIn(BAM_PLATFORMS, default="ILLUMINA"),
    # Offset for quality scores in FASTQ files.
    "QualityOffset": ValueIn((33, 64, "Solexa"), default=33),
    # Split a lane into multiple entries, one for each (pair of) file(s)
    "AdapterRemoval": {
        "--adapter1": IsStr,
        "--adapter2": IsStr,
        "--adapter-list": IsStr,
        "--maxns": IsUnsignedInt,
        "--minquality": IsUnsignedInt,
        "--trimns": IsNone,
        "--trimqualities": IsNone,
        "--collapse": IsNone,
        "--mm": Or(IsFloat, IsUnsignedInt),
        "--minlength": IsUnsignedInt,
        "--maxlength": IsUnsignedInt,
        "--minalignmentlength": IsUnsignedInt,
        "--minadapteroverlap": IsUnsignedInt,
        "--shift": IsUnsignedInt,
        "--qualitymax": IsUnsignedInt,
        "--mate-separator": IsStr,
        "--trimwindows": Or(IsInt, IsFloat),
        "--preserve5p": IsNone,
        "--collapse-deterministic": IsNone,
        "--collapse-conservatively": IsNone,
    },
    # Which aliger/mapper to use (BWA/Bowtie2)
    "Aligners": {
        "Program": ValueIn(("BWA", "Bowtie2"), default="BWA"),
        "BWA": {
            # Mapping algorithm; availability depends on BWA version
            "Algorithm": StringIn(("backtrack", "mem", "bwasw"), default="backtrack"),
            # Minimum mapping quality (PHREAD) of reads to retain
            "MinQuality": IsUnsignedInt(default=0),
            # Remove unmapped reads or not
            "FilterUnmappedReads": IsBoolean(default=True),
            # Use seed region during mapping
            # Verbose name for command-line option "-l 65535"
            "UseSeed": IsBoolean(default=True),
            # Any number of user specific options
            StringStartsWith("-"): Or(
                IsListOf(IsStr, IsInt, IsFloat), Or(IsStr, IsInt, IsFloat, IsNone)
            ),
        },
        "Bowtie2": {
            # Minimum mapping quality (PHREAD) of reads to retain
            "MinQuality": IsUnsignedInt(default=0),
            # Remove unmapped reads or not
            "FilterUnmappedReads": IsBoolean(default=True),
            # Any number of user specific options
            StringStartsWith("-"): Or(
                IsListOf(IsStr, IsInt, IsFloat), Or(IsStr, IsInt, IsFloat, IsNone)
            ),
        },
    },
    "mapDamage": {
        # Tabulation options
        "--downsample": Or(IsUnsignedInt, IsFloat),
        "--length": IsUnsignedInt,
        "--around": IsUnsignedInt,
        "--min-basequal": IsUnsignedInt,
        # Plotting options
        "--ymax": IsFloat,
        "--readplot": IsUnsignedInt,
        "--refplot": IsUnsignedInt,
        # Model options
        "--rand": IsUnsignedInt,
        "--burn": IsUnsignedInt,
        "--adjust": IsUnsignedInt,
        "--iter": IsUnsignedInt,
        "--forward": IsNone,
        "--reverse": IsNone,
        "--var-disp": IsNone,
        "--jukes-cantor": IsNone,
        "--diff-hangs": IsNone,
        "--fix-nicks": IsNone,
        "--use-raw-nick-freq": IsNone,
        "--single-stranded": IsNone,
        "--seq-length": IsUnsignedInt,
    },
    # Exclude READ_TYPES from alignment/analysis
    "ExcludeReads": _VALID_EXCLUDE_DICT,
    # Features of pipeline
    "Features": _VALID_FEATURES_DICT,
}


_VALIDATION = {
    "Options": _VALIDATION_OPTIONS,
    "Genomes": {
        _VALID_GENOME_NAME: {
            "Path": _VALID_GENOME_PATH,
        },
    },
    _VALID_FILENAME: {  # Group
        _VALID_FILENAME: {  # Sample
            _VALID_FILENAME: {  # Library
                _VALID_FILENAME: Or(IsStr, IsDictOf(IsStr, IsStr)),
                "Options": WithoutDefaults(_VALIDATION_OPTIONS),
            },
            "Options": WithoutDefaults(_VALIDATION_OPTIONS),
        },
        "Options": WithoutDefaults(_VALIDATION_OPTIONS),
    },
}


def _mangle_makefile(makefile, pipeline_variant):
    options = makefile.pop("Options")
    genomes = makefile.pop("Genomes")
    makefile = {
        "Options": options,
        "Genomes": genomes,
        "Samples": _flatten_groups(options, makefile),
    }

    _mangle_genomes(makefile, require_genome=pipeline_variant != "trim")
    _mangle_lanes(makefile)
    _mangle_tags(makefile)

    _split_lanes_by_filenames(makefile)

    # FIXME: Remove from rest of code-base
    makefile["Targets"] = {
        sample: {sample: libraries} for sample, libraries in makefile["Samples"].items()
    }

    # FIXME: Update from rest of code-base
    makefile["Prefixes"] = makefile["Genomes"]

    return makefile


def _flatten_groups(options, groups):
    results = {}
    for group, samples in groups.items():
        _merge_options(options, group, samples)

        for name, sample in samples.items():
            if name in results:
                raise MakefileError("Multiple samples with name {!r}".format(name))

            results[name] = sample

    return results


def _merge_options(options, group, samples):
    def _combine_options(options, items, path, valid_features=_VALID_FEATURES_DICT):
        item_options = items.pop("Options", {})
        if not item_options:
            return options

        features = item_options.get("Features", {})
        invalid_features = features.keys() - valid_features
        if invalid_features:
            raise MakefileError(
                "Cannot enable/disable %s on a per-library basis at %s"
                % (", ".join(map(repr, invalid_features)), " :: ".join(path))
            )

        # Fill out missing values using those of prior levels
        return fill_dict(destination=item_options, source=options)

    # Options common to a specific group
    group_options = _combine_options(options, samples, (group,))
    for sample, libraries in samples.items():
        # Options common to a specific sample in a group
        sample_options = _combine_options(group_options, libraries, (group, sample))

        for library, lanes in libraries.items():
            # Options common to a specific library in a sample
            lanes["Options"] = _combine_options(
                sample_options,
                lanes,
                (group, sample, library),
                _VALID_LIBRARY_LEVEL_FEATURES,
            )


def _mangle_genomes(makefile, require_genome):
    genomes = makefile["Genomes"]
    for name in tuple(genomes):
        if "*" in name[:-1]:
            raise MakefileError(
                "The character '*' is not allowed in Genome names; if you wish to "
                "select multiple .fasta files using a search-string, then use the "
                "genome name '%s*' instead and specify the wildcards in the 'Path'."
                % (name.replace("*", ""))
            )
        elif name.endswith("*"):
            filename = genomes.pop(name)["Path"]
            for name, filename in _glob_genomes(filename):
                if name in genomes:
                    raise MakefileError("Multiple genomes with the name %s" % name)

                genomes[name] = {
                    "Name": name,
                    "Path": filename,
                }
        else:
            genomes[name]["Name"] = name

    if require_genome and not genomes:
        raise MakefileError("At least one genome must be specified")


def _glob_genomes(pattern):
    filename = None
    for filename in glob.iglob(pattern):
        name = os.path.basename(filename).split(".")[0]
        _VALID_GENOME_NAME(("Genomes", name), name)

        yield (name, filename)

    if filename is None:
        raise MakefileError("Did not find any genomes using wildcards %r" % (pattern,))


def _mangle_lanes(makefile):
    formatter = string.Formatter()
    for (sample_name, libraries) in makefile["Samples"].items():
        for (library_name, lanes) in libraries.items():
            options = lanes.pop("Options")

            for (lane, data) in lanes.items():
                path = (sample_name, library_name, lane)

                _validate_lane_paths(data, path, formatter)

                lane_type = _determine_lane_type(data, path)
                if lane_type == "Trimmed" and options["QualityOffset"] != 33:
                    raise MakefileError(
                        "Pre-trimmed data must have quality offset 33 (Phred+33). "
                        "Please convert your FASTQ files using e.g. seqtk before "
                        "continuing: {}".format(" :: ".join(path))
                    )

                lanes[lane] = {"Type": lane_type, "Data": data, "Options": options}


def _validate_lane_paths(data, path, fmt):
    filenames = []
    if isinstance(data, str):
        filenames.append(data)
    elif isinstance(data, dict):
        filenames.extend(iter(data.values()))

    for filename in filenames:
        try:
            fields = tuple(fmt.parse(filename))
        except ValueError as error:
            raise MakefileError(
                "Error parsing path specified at %r; %s; note "
                "that the characters '}' and '{' should only "
                "be used as part of the key '{Pair}', in "
                "order to specify the mate identifier: %r"
                % (" :: ".join(path), error, filename)
            )

        for _, key, _, _ in fields:
            if key not in (None, "Pair"):
                raise MakefileError(
                    "Invalid path specified at %r; only the "
                    "key '{Pair}' is allowed, to specify the "
                    "mate 1 / 2 identifier, but the key "
                    "'{%s}' was found in the path: %r"
                    % (" :: ".join(path), key, filename)
                )


def _determine_lane_type(data, path):
    if isinstance(data, str):
        return "Raw"

    bad_keys = data.keys() - _READ_TYPES
    if bad_keys:
        raise MakefileError(
            "Error at %s: Keys must be any of Paired, Single, Collapsed, "
            "CollapsedTruncated, or Singleton, but found %s"
            % (" :: ".join(path), ", ".join(map(repr, bad_keys))),
        )

    for (key, files) in data.items():
        is_paired = paths.is_paired_end(files)

        if is_paired and key != "Paired":
            raise MakefileError(
                "Error at %s: Path includes {Pair} key, but read-type "
                "is not Paired: %r" % (" :: ".join(path + (key,)), files)
            )
        elif not is_paired and key == "Paired":
            raise MakefileError(
                "Error at %s: Paired pre-trimmed reads specified, but path "
                "does not contain {Pair} key: %r" % (" :: ".join(path + (key,)), files)
            )

    return "Trimmed"


def _mangle_tags(makefile):
    for (sample, libraries) in makefile["Samples"].items():
        for (library, barcodes) in libraries.items():
            for (barcode, record) in barcodes.items():
                record["Tags"] = {
                    # FIXME: Remove
                    "Target": sample,
                    "ID": library,
                    "SM": sample,
                    "LB": library,
                    "PU": barcode,
                    "DS": "NA",
                    "Folder": "NA",
                    "PG": record["Options"]["Aligners"]["Program"],
                    "PL": record["Options"]["Platform"].upper(),
                }


def _split_lanes_by_filenames(makefile):
    for (sample, library, barcode, record) in _iterate_over_records(makefile):
        if record["Type"] == "Raw":
            path = (sample, library, barcode)
            filenames = paths.collect_files(path, record["Data"])
            filename_keys = sorted(filenames)  # Either ["SE"] or ["PE_1", "PE_2"]

            library = makefile["Samples"][sample][library]
            template = library.pop(barcode)

            input_files = [filenames[key] for key in filename_keys]
            input_files_iter = itertools.zip_longest(*input_files)
            for (index, filenames) in enumerate(input_files_iter, start=1):
                current = copy.deepcopy(template)

                assert len(filenames) == len(filename_keys), filenames
                current["Data"] = dict(zip(filename_keys, filenames))

                # Save a summary of file paths in description (DS) tag
                current["Tags"]["DS"] = _summarize_filenames(
                    filenames, show_differences=True
                )

                # Save a summary of file names for use as temporary folders
                current["Tags"]["Folder"] = _summarize_filenames(
                    [os.path.basename(filename) for filename in filenames]
                ).replace("?", "x")

                # ':' is disallowed in barcodes and so are safe for generated names
                new_barcode = "%s::%03i" % (barcode, index)
                assert new_barcode not in library, (new_barcode, library)
                library[new_barcode] = current


def _summarize_filenames(filenames, show_differences=False):
    combined_filenames = get_files_glob(filenames, show_differences=show_differences)
    if combined_filenames is None:
        return ";".join(filenames)

    return combined_filenames


def _validate_makefiles(makefiles):
    for makefile in makefiles:
        _validate_makefile_adapters(makefile)
    _validate_makefiles_duplicate_samples(makefiles)
    _validate_makefiles_duplicate_files(makefiles)
    _validate_prefixes(makefiles)

    return makefiles


def _validate_makefile_adapters(makefile):
    """Checks for the default adapter sequences specified in the wrong
    orientation for AdapterRemoval, which is a typical mistake when using
    the --pcr2 option.
    """
    # The non-reverse complemented mate 2 adapter, as seen in raw FASTQ reads
    adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

    tests = {
        # --pcr2 expects the reverse complement of the mate 2 adapter seq.
        "--pcr2": adapter_2,
        # --adapter2 (AdapterRemoval v2) expects the regular sequence
        "--adapter2": sequences.reverse_complement(adapter_2),
    }

    def check_options(options, results):
        for key, value in tests.items():
            if options.get(key) == value:
                results[key] = True

    results = dict.fromkeys(tests, False)
    for (_, _, _, record) in _iterate_over_records(makefile):
        adapterrm_opt = record.get("Options", {}).get("AdapterRemoval", {})
        check_options(adapterrm_opt, results)

    adapterrm_opt = makefile.get("Options", {}).get("AdapterRemoval", {})
    check_options(adapterrm_opt, results)

    if any(results.values()):
        logger = logging.getLogger(__name__)
        logger.warn(
            "An adapter specified for AdapterRemoval corresponds to the default "
            "sequence, but is reverse complemented. Please make sure that this is "
            "intended! "
        )

        if results["--pcr2"]:
            logger.warn(
                "For --pcr2, the sequence given should be the "
                "reverse complement of the sequence observed in the "
                "mate 2 FASTQ file."
            )

        if results["--adapter2"]:
            logger.warn(
                "For --adapter2 (AdapterRemoval v2, only) the value "
                "should be exactly as observed in the FASTQ reads."
            )


def _validate_makefiles_duplicate_files(makefiles):
    filenames = collections.defaultdict(list)
    for makefile in makefiles:
        iterator = _iterate_over_records(makefile)
        for (sample, library, barcode, record) in iterator:
            for realpath in map(os.path.realpath, record["Data"].values()):
                filenames[realpath].append((sample, library, barcode))

    has_overlap = {}
    for (filename, records) in filenames.items():
        if len(records) > 1:
            has_overlap[filename] = list(set(records))

    logger = logging.getLogger(__name__)
    by_records = sorted(zip(list(has_overlap.values()), list(has_overlap.keys())))
    for (records, pairs) in itertools.groupby(by_records, lambda x: x[0]):
        pairs = list(pairs)
        description = _describe_files_in_multiple_records(records, pairs)

        if len(set(record[0] for record in records)) != len(records):
            message = "Path included multiple times in sample:\n"
            raise MakefileError(message + description)
        else:
            logger.warn("WARNING: Path included in multiple samples:\n%s", description)


def _describe_files_in_multiple_records(records, pairs):
    descriptions = []
    for (index, record) in enumerate(sorted(records), start=1):
        descriptions.append(
            "\t- Record {0}: Name: {1},  Sample: {2},  "
            "Library: {3},  Barcode: {4}".format(index, *record)
        )

    for (index, (_, filename)) in enumerate(sorted(pairs), start=1):
        message = "\t- Canonical path {0}: {1}"
        descriptions.append(message.format(index, filename))

    return "\n".join(descriptions)


def _validate_makefiles_duplicate_samples(makefiles):
    samples = set()
    for makefile in makefiles:
        for sample in samples - makefile["Samples"].keys():
            raise MakefileError(
                "Sample name '%s' used multiple times; output files would be clobbered!"
                % (sample,)
            )
        samples.update(makefile["Samples"])


def _validate_prefixes(makefiles):
    logger = logging.getLogger(__name__)
    already_validated = {}
    logger.info("Validating FASTA files")
    for makefile in makefiles:
        for prefix in makefile["Prefixes"].values():
            path = prefix["Path"]
            if path in already_validated:
                prefix["IndexFormat"] = already_validated[path]["IndexFormat"]
                continue

            # Must be set to a valid value, even if FASTA file does not exist
            prefix["IndexFormat"] = ".bai"

            if not os.path.exists(path):
                logger.error("Reference FASTA file does not exist: %r", path)
                continue
            elif not os.path.exists(path + ".fai"):
                logger.info("Indexing FASTA at %r", path)

            try:
                contigs = FASTA.index_and_collect_contigs(path)
            except FASTAError as error:
                raise MakefileError("Error indexing FASTA:\n %s" % (error,))

            if max(contigs.values()) > _BAM_MAX_SEQUENCE_LENGTH:
                logger.warn(
                    "FASTA file %r contains sequences longer "
                    "than %i! CSI index files will be used instead "
                    "of BAI index files.",
                    path,
                    _BAM_MAX_SEQUENCE_LENGTH,
                )
                prefix["IndexFormat"] = ".csi"

            already_validated[path] = prefix


def _iterate_over_records(makefile):
    for (sample, libraries) in tuple(makefile["Samples"].items()):
        for (library, barcodes) in tuple(libraries.items()):
            for (barcode, record) in tuple(barcodes.items()):
                yield sample, library, barcode, record
