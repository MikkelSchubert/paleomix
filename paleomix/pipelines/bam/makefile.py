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
from typing import Any, Dict, Iterable, Tuple

import paleomix.common.sequences as sequences
import paleomix.pipelines.bam.paths as paths
from paleomix.common.bamfiles import BAM_PLATFORMS
from paleomix.common.fileutils import get_files_glob
from paleomix.common.formats.fasta import FASTA
from paleomix.common.makefile import (
    REQUIRED_VALUE,
    And,
    IsAny,
    IsBoolean,
    IsFloat,
    IsInt,
    IsListOf,
    IsNone,
    IsStr,
    IsUnsignedInt,
    MakefileError,
    Not,
    Or,
    StringIn,
    StringStartsWith,
    ValueIn,
    ValuesIntersect,
    ValuesSubsetOf,
    WithoutDefaults,
    process_makefile,
    read_makefile,
)
from paleomix.common.utilities import fill_dict

_READ_TYPES = set(("Single", "Singleton", "Collapsed", "CollapsedTruncated", "Paired"))

# The maximum reference sequence length supported by the BAI index format:
#   https://samtools.github.io/hts-specs/SAMv1.pdf
_BAM_MAX_SEQUENCE_LENGTH = 2 ** 29 - 1


def read_makefiles(filenames):
    logger = logging.getLogger(__name__)

    makefiles = []
    for filename in filenames:
        logger.info("Reading makefile %r", filename)

        data = read_makefile(filename, MAKEFILE_SPECIFICATION)
        options = data.pop("Options")
        genomes = data.pop("Genomes")

        makefile = {
            "Filename": filename,
            "Options": options,
            "Genomes": _postprocess_genomes(genomes),
            "Samples": _postprocess_samples(data, options),
        }

        # FIXME: Remove from rest of code-base
        makefile["Targets"] = {
            sample: {sample: libraries}
            for sample, libraries in makefile["Samples"].items()
        }

        # FIXME: Remove from rest of code-base
        makefile["Prefixes"] = makefile["Genomes"]

        makefiles.append(makefile)

    return validate_makefiles(makefiles)


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
    Not(StringIn(["Options"])),
)

# Valid paths for genomes; avoids some problems with e.g. Bowtie2
_VALID_GENOME_PATH = And(
    IsStr(), Not(ValuesIntersect('\\:?"<>|() \t\n\v\f\r')), default=REQUIRED_VALUE
)

# Valid strings for samples / libraries / lanes
_VALID_FILENAME = _alphanum_check(whitelist="._-", min_len=2)


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
    "ExcludeReads": {
        "Single": IsBoolean(default=False),
        "Collapsed": IsBoolean(default=False),
        "CollapsedTruncated": IsBoolean(default=False),
        "Paired": IsBoolean(default=False),
        "Singleton": IsBoolean(default=False),
    },
    # Features of pipeline
    "Features": {
        "Coverage": IsBoolean(default=True),
        "Depths": IsBoolean(default=True),
        "mapDamage": StringIn(("rescale", "model", "plot", True, False), default=False),
        "PCRDuplicates": StringIn((True, False, "mark", "filter"), default="filter"),
        "Summary": IsBoolean(default=True),
    },
}

# validation of a complex lane, containing trimmed reads and/or options
_VALIDATE_LANE = {
    "Single": IsStr,
    "Collapsed": IsStr,
    "CollapsedTruncated": IsStr,
    "Paired": IsStr,
    "Singleton": IsStr,
    "Raw": IsStr,
    "Options": WithoutDefaults(_VALIDATION_OPTIONS),
}

MAKEFILE_SPECIFICATION = {
    "Options": _VALIDATION_OPTIONS,
    "Genomes": {
        _VALID_GENOME_NAME: {
            "Path": _VALID_GENOME_PATH,
        },
    },
    _VALID_FILENAME: {  # Group
        _VALID_FILENAME: {  # Sample
            _VALID_FILENAME: {  # Library
                # Validation of lanes is performed in `_postprocess_samples`
                _VALID_FILENAME: IsAny,
                "Options": WithoutDefaults(_VALIDATION_OPTIONS),
            },
            "Options": WithoutDefaults(_VALIDATION_OPTIONS),
        },
        "Options": WithoutDefaults(_VALIDATION_OPTIONS),
    },
}


########################################################################################
# Post processing of user defined target genomes


def _postprocess_genomes(genomes):
    result = {}
    for name, values in genomes.items():
        if "*" in name[:-1]:
            raise MakefileError(
                "The character '*' is not allowed in Genome names; if you wish to "
                "select multiple .fasta files using a search-string, then use the "
                "genome name '%s*' instead and specify the wildcards in the 'Path'."
                % (name.replace("*", ""))
            )
        elif name.endswith("*"):
            for name, filename in _glob_genomes(values["Path"]):
                if name in result:
                    raise MakefileError(f"Multiple genomes named {name}")

                result[name] = {
                    "Name": name,
                    "Path": filename,
                }
        elif name in result:
            raise MakefileError(f"Multiple genomes named {name}")
        else:
            result[name] = {
                "Name": name,
                "Path": values["Path"],
            }

    return result


def _glob_genomes(pattern):
    filename = None
    for filename in glob.iglob(pattern):
        name = os.path.basename(filename).split(".")[0]
        _VALID_GENOME_NAME(("Genomes", name), name)

        yield (name, filename)

    if filename is None:
        raise MakefileError(f"Did not find any genomes using wildcards {pattern!r}")


########################################################################################
# Post processing of user samples / sample groups


def _postprocess_samples(data, global_options):
    top_level_features = ("Coverage", "Depths", "mapDamage", "PCRDuplicates", "Summary")

    for (group, samples) in tuple(data.items()):
        # Options common to a specific group
        group_options = _combine_options(
            options=global_options,
            data=samples,
            path=(group,),
            valid_features=top_level_features,
        )

        for sample, libraries in samples.items():
            # Options common to a specific sample in a group
            sample_options = _combine_options(
                options=group_options,
                data=libraries,
                path=(group, sample),
                valid_features=top_level_features,
            )

            for library, lanes in libraries.items():
                # Options common to a specific library in a sample
                library_options = _combine_options(
                    options=sample_options,
                    data=lanes,
                    path=(group, sample, library),
                    valid_features=("mapDamage", "PCRDuplicates"),
                )

                for barcode, record in lanes.items():
                    path = (group, sample, library, barcode)

                    # lane data must be normalized since there are 2 possible forms
                    record = _normalize_lane(record, path)

                    # once the lanes are normalized, we generate lane-specific options
                    record["Options"] = _combine_options(
                        options=library_options,
                        data=record,
                        path=path,
                        valid_features=(),
                    )

                    # Record information used for tagging/naming files
                    _add_lane_tags(record, group, sample, library, barcode)

                    lanes[barcode] = record

    # Flatten groups and split lanes by filenames
    return _finalize_samples(data)


def _combine_options(
    options: Dict[str, Any],
    data: Dict[str, Any],
    path: Tuple[str, ...],
    valid_features: Iterable[str],
) -> Dict[str, Any]:
    item_options = data.pop("Options", {})
    if not item_options:
        return options

    features = item_options.get("Features", {})
    invalid_features = features.keys() - valid_features
    if invalid_features:
        raise MakefileError(
            "Cannot override %s at %s"
            % (", ".join(map(repr, invalid_features)), _path_to_str(path))
        )

    # Fill out missing values using those of prior levels
    return fill_dict(destination=item_options, source=options)


def _normalize_lane(data, path) -> Dict[str, Any]:
    if isinstance(data, str):
        # The simple case: raw data without options
        return {"Type": "Raw", "Data": data}

    # the structure needs to be validated here, since the specification uses an IsAny
    data = process_makefile(
        data=data,
        specification=_VALIDATE_LANE,
        path=path,
        apply_defaults=False,
    )

    options = data.pop("Options", None)

    if "Raw" in data:
        # it doesn't make sense to have both trimmed and untrimmed reads in one lane
        if data.keys() & _READ_TYPES:
            raise MakefileError(f"both raw and trimmed reads at {_path_to_str(path)}")

        ((type_, data),) = data.items()
    else:
        type_ = "Trimmed"

    data = {"Type": type_, "Data": data}
    if options is not None:
        data["Options"] = options

    return data


def _add_lane_tags(record, group, sample, library, barcode):
    record["Tags"] = {
        "Target": sample,  # FIXME: Remove
        "Group": group,
        "ID": library,
        "SM": sample,
        "LB": library,
        "PU": barcode,
        "DS": "NA",
        "Folder": "NA",
        "PG": record["Options"]["Aligners"]["Program"],
        "PL": record["Options"]["Platform"].upper(),
    }


def _finalize_samples(data):
    results = {}
    for group, samples in data.items():
        for sample, libraries in samples.items():
            if sample in results:
                raise MakefileError(f"Multiple samples named {sample!r}")

            split_libraries = {}
            for library, lanes in libraries.items():
                split_lanes = {}
                for barcode, record in lanes.items():
                    split_lanes.update(
                        _split_lanes_by_filenames(
                            group=group,
                            sample=sample,
                            library=library,
                            barcode=barcode,
                            record=record,
                        )
                    )

                split_libraries[library] = split_lanes
            results[sample] = split_libraries

    return results


def _split_lanes_by_filenames(group, sample, library, barcode, record):
    if record["Type"] == "Raw":
        path = (group, sample, library, barcode)
        filenames = paths.collect_files(path, record["Data"])
        filename_keys = sorted(filenames)  # Either ["SE"] or ["PE_1", "PE_2"]

        input_files = [filenames[key] for key in filename_keys]
        input_files_iter = itertools.zip_longest(*input_files)
        for (index, filenames) in enumerate(input_files_iter, start=1):
            current = copy.deepcopy(record)

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
            yield (f"{barcode}::{index:03}", current)
    else:
        yield (barcode, record)


def _summarize_filenames(filenames, show_differences=False):
    combined_filenames = get_files_glob(filenames, show_differences=show_differences)
    if combined_filenames is None:
        return ";".join(filenames)

    return combined_filenames


########################################################################################


def validate_makefiles(makefiles):
    for makefile in makefiles:
        _validate_makefile_options(makefile)
    _validate_makefiles_duplicate_samples(makefiles)
    _validate_makefiles_duplicate_files(makefiles)
    _validate_prefixes(makefiles)

    return makefiles


def _validate_makefile_options(makefile):
    for (sample, library, barcode, record) in _iterate_over_records(makefile):
        path = (record["Tags"]["Group"], sample, library, barcode)

        if record["Type"] == "Trimmed":
            if record["Options"]["QualityOffset"] != 33:
                raise MakefileError(
                    "Pre-trimmed data must have quality offset 33 (Phred+33). "
                    "Please convert your FASTQ files using e.g. seqtk before "
                    "continuing: {}".format(_path_to_str(path))
                )

            for (key, files) in record["Data"].items():
                is_paired = paths.is_paired_end(files)

                if is_paired and key != "Paired":
                    raise MakefileError(
                        "Error at %s: Path includes {Pair} key, but read-type "
                        "is not Paired: %r" % (_path_to_str(path + (key,)), files)
                    )
                elif not is_paired and key == "Paired":
                    raise MakefileError(
                        "Error at %s: Paired pre-trimmed reads specified, but path "
                        "does not contain {Pair} key: %r"
                        % (_path_to_str(path + (key,)), files)
                    )

        _validate_makefile_adapters(record, path)


def _validate_makefile_adapters(record, path):
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

    results = {}
    options = record["Options"]
    for key, value in tests.items():
        results[key] = options.get(key, "").upper() == value

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
        for (sample, library, _, record) in _iterate_over_records(makefile):
            # use user defined barcode rather than auto-numbered barcode
            barcode = record["Tags"]["PU"]

            for realpath in map(os.path.realpath, record["Data"].values()):
                filenames[realpath].append((sample, library, barcode))

    has_overlap = {}
    for (filename, records) in filenames.items():
        if len(records) > 1:
            has_overlap[filename] = list(set(records))

    logger = logging.getLogger(__name__)
    by_records = sorted(zip(list(has_overlap.values()), list(has_overlap.keys())))
    for (records, pairs) in itertools.groupby(by_records, lambda x: x[0]):
        description = _describe_files_in_multiple_records(records, pairs)

        if len(set(record[0] for record in records)) != len(records):
            message = "FASTQ files are used multiple times in sample:\n"
            raise MakefileError(message + description)
        else:
            logger.warn("WARNING: Path included in multiple samples:\n%s", description)


def _describe_files_in_multiple_records(records, pairs):
    lines = []
    prefix = "Filename"
    for (_, filename) in sorted(pairs):
        lines.append("  {0} {1}".format(prefix, filename))
        prefix = " " * len(prefix)

    prefix = "Found at"
    for record in sorted(records):
        # FIXME: Show the glob that found the above files
        lines.append("  {0} {1} :: {2} :: {3}".format(prefix, *record))
        prefix = " " * len(prefix)

    return "\n".join(lines)


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
            except Exception as error:
                raise MakefileError("Error reading/indexing FASTA: %s" % (error,))

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


def _path_to_str(path):
    return " :: ".join(path)


def _iterate_over_records(makefile):
    for (sample, libraries) in tuple(makefile["Samples"].items()):
        for (library, barcodes) in tuple(libraries.items()):
            for (barcode, record) in tuple(barcodes.items()):
                yield sample, library, barcode, record
