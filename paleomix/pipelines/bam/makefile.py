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
from __future__ import annotations

import collections
import glob
import itertools
import logging
import os
import re
import string
from typing import Any, Dict, Iterable, Optional, Tuple

from paleomix.common import sequences
from paleomix.common.bamfiles import BAM_PLATFORMS
from paleomix.common.fileutils import get_files_glob
from paleomix.common.formats.fasta import FASTA
from paleomix.common.makefile import (
    REQUIRED_VALUE,
    And,
    DeprecatedOption,
    FASTQPath,
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
    SpecTree,
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
_BAM_MAX_SEQUENCE_LENGTH = 2**29 - 1


def read_makefiles(filenames: Iterable[str], pipeline_variant: str = "bam"):
    logger = logging.getLogger(__name__)

    makefiles = []
    for filename in filenames:
        logger.info("Reading makefile %r", filename)

        data = read_makefile(filename, MAKEFILE_SPECIFICATION)
        if not isinstance(data, dict):
            raise AssertionError("makefile is not a dict")

        options = data.pop("Options")
        genomes = data.pop("Genomes")

        makefiles.append(
            {
                "Filename": filename,
                "Options": options,
                "Genomes": _postprocess_genomes(genomes),
                "Samples": _postprocess_samples(data, options),
            }
        )

    return finalize_makefiles(makefiles, pipeline_variant)


def _alphanum_check(whitelist, min_len=1):
    description = "characters a-z, A-Z, 0-9%s allowed"
    description %= (", and %r" % whitelist,) if whitelist else ""

    whitelist += string.ascii_letters + string.digits

    return And(
        IsStr(min_len=min_len), ValuesSubsetOf(whitelist, description=description)
    )


# Valid names for genomes
_VALID_GENOME_NAME: SpecTree = And(
    _alphanum_check(whitelist="._-*"),
    Not(ValueIn(["Options"])),
)

# Valid paths for genomes; avoids some problems with e.g. Bowtie2
_VALID_GENOME_PATH: SpecTree = And(
    IsStr(), Not(ValuesIntersect('\\:?"<>|() \t\n\v\f\r')), default=REQUIRED_VALUE
)

# Valid strings for samples / libraries / lanes
_VALID_FILENAME: SpecTree = _alphanum_check(whitelist="._-", min_len=2)


_VALIDATION_OPTIONS: SpecTree = {
    # Sequencing platform, used to tag read-groups.
    "Platform": ValueIn(BAM_PLATFORMS, default="ILLUMINA"),
    # Offset for quality scores in FASTQ files.
    "QualityOffset": ValueIn((33, 64, "Solexa"), default=33),
    # Split a lane into multiple entries, one for each (pair of) file(s)
    "AdapterRemoval": {
        "Version": ValueIn((2, 3), default=2),
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
        "--trim5p": Or(IsInt, IsListOf(IsInt)),
        "--trim3p": Or(IsInt, IsListOf(IsInt)),
        StringStartsWith("--"): Or(IsStr, IsInt, IsFloat, IsNone),
    },
    # Which aliger/mapper to use (BWA/Bowtie2)
    "Aligners": {
        "Program": ValueIn(("BWA", "Bowtie2"), default="BWA"),
        "BWA": {
            # Mapping algorithm; availability depends on BWA version
            "Algorithm": ValueIn(
                ("backtrack", "mem", "mem2", "bwasw"),
                default="mem",
            ),
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
        "mapDamage": ValueIn(("rescale", "model", "plot", True, False), default=False),
        "PCRDuplicates": ValueIn((True, False, "mark", "filter"), default="filter"),
        # TODO: Statistics to be combined into new report (HTML + JSON?)
        "Coverage": DeprecatedOption(IsBoolean(default=True)),
        "Depths": DeprecatedOption(IsBoolean(default=True)),
        "Summary": DeprecatedOption(IsBoolean(default=True)),
    },
}

# validation of a complex lane, containing trimmed reads and/or options
_VALIDATE_LANE: SpecTree = {
    "Single": FASTQPath,
    "Collapsed": FASTQPath,
    "CollapsedTruncated": FASTQPath,
    "Paired": FASTQPath(paired_end=True),
    "Singleton": FASTQPath,
    "Untrimmed": FASTQPath(paired_end=None),
    "Options": WithoutDefaults(_VALIDATION_OPTIONS),
}

MAKEFILE_SPECIFICATION: SpecTree = {
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

    for group, samples in tuple(data.items()):
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

                    # Split a trimmed/untrimmed lane into one record per input file
                    lanes[barcode] = _split_lane(record, path, library_options)

    return data


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

    if "mapDamage" in item_options and "mapDamage" not in valid_features:
        raise MakefileError(f"Cannot set mapDamage options at {_path_to_str(path)}")

    # Fill out missing values using those of prior levels
    return fill_dict(destination=item_options, source=options)


def _split_lane(data, path, options):
    if isinstance(data, str):
        data = {"Untrimmed": data}

    # the structure needs to be validated here, since the specification uses an IsAny
    data = process_makefile(
        data=data,
        specification=_VALIDATE_LANE,
        path=path,
        apply_defaults=False,
    )

    # Generate final options for this lane
    options = _combine_options(
        options=options,
        data=data,
        path=path,
        valid_features=(),
    )

    return [
        {
            "Path": files,
            "Type": read_type,
            "Shortname": _filenames_to_shortname(files),
            "Options": options,
        }
        for read_type, files in _collect_files_and_split_lane(data, path)
    ]


def _collect_files_and_split_lane(data, path):
    if "Untrimmed" in data and len(data) > 1:
        raise MakefileError(f"both untrimmed and trimmed reads at {_path_to_str(path)}")

    for read_type, filepath in data.items():
        if read_type == "Untrimmed":
            for files in _collect_files(path, filepath):
                yield read_type, files
        elif read_type == "Paired":
            yield (
                read_type,
                (
                    FASTQPath.format(filepath, 1),
                    FASTQPath.format(filepath, 2),
                ),
            )
        else:
            yield read_type, (filepath, None)


def _filenames_to_shortname(filenames):
    if not (1 <= len(filenames) <= 2):
        raise ValueError(filenames)

    basenames = []
    for filename in filenames:
        if filename is not None:
            basenames.append(os.path.basename(filename))

    filename = get_files_glob(basenames)
    if filename is None:
        raise ValueError(filenames)

    return filename.replace("?", "x")


def _collect_files(path, template) -> Iterable[Tuple[str, Optional[str]]]:
    if FASTQPath.is_paired_end(template):
        files_1 = _sorted_glob(FASTQPath.format(template, 1))
        files_2 = _sorted_glob(FASTQPath.format(template, 2))

        if len(files_1) != len(files_2):
            raise MakefileError(
                "Unequal number of mate 1 and mate 2 files found at %r; found %i "
                "mate 1 files and %i mate 2 files; specified in makefile at %r. "
                % (template, len(files_1), len(files_2), _path_to_str(path))
            )
        elif not (files_1 and files_2):
            return [(template, None)]

        return zip(files_1, files_2)
    else:
        files = _sorted_glob(template)
        if not files:
            return [(template, None)]

        return [(filename, None) for filename in files]


def _sorted_glob(filename):
    if _GLOB_MAGIC.search(filename):
        return sorted(glob.iglob(filename))

    return [filename]


# based on check in `glob`
_GLOB_MAGIC = re.compile("[*?[]")

########################################################################################


def finalize_makefiles(makefiles, pipeline_variant):
    sample_names = set()
    duplicate_samples = set()
    for makefile in makefiles:
        results = {}
        # Groups are used to structure the YAML file and can be discarded for simplicity
        for samples in makefile["Samples"].values():
            for sample, libraries in samples.items():
                if sample in sample_names:
                    duplicate_samples.add(sample)

                sample_names.add(sample)
                results[sample] = libraries

        makefile["Samples"] = results

    if duplicate_samples:
        log = logging.getLogger(__name__)
        log.error("One or more sample names have been used multiple times:")
        for idx, sample in enumerate(sorted(duplicate_samples), start=1):
            log.error("  %i. %s", idx, sample)
        log.error("All samples must have a unique name")

        raise MakefileError("Duplicate sample names found")

    for makefile in makefiles:
        _validate_makefile_options(makefile)
    _validate_makefiles_duplicate_files(makefiles)
    _validate_prefixes(makefiles, pipeline_variant)

    return makefiles


def _validate_makefile_options(makefile):
    for sample, library, barcode, record in _iterate_over_records(makefile):
        path = (sample, library, barcode)

        if record["Type"] != "Untrimmed":
            if record["Options"]["QualityOffset"] != 33:
                raise MakefileError(
                    "Pre-trimmed data must have quality offset 33 (Phred+33). "
                    "Please convert your FASTQ files using e.g. seqtk before "
                    "continuing: {}".format(_path_to_str(path))
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
        for sample, library, barcode, record in _iterate_over_records(makefile):
            for filepath in record["Path"]:
                if filepath is not None:
                    realpath = os.path.realpath(filepath)

                    filenames[realpath].append((sample, library, barcode))

    has_overlap = {}
    for filename, records in filenames.items():
        if len(records) > 1:
            has_overlap[filename] = list(set(records))

    logger = logging.getLogger(__name__)
    by_records = sorted(zip(list(has_overlap.values()), list(has_overlap.keys())))
    for records, pairs in itertools.groupby(by_records, lambda x: x[0]):
        description = _describe_files_in_multiple_records(records, pairs)

        if len(set(record[0] for record in records)) != len(records):
            message = "FASTQ files are used multiple times in sample:\n"
            raise MakefileError(message + description)
        else:
            logger.warn("WARNING: Path included in multiple samples:\n%s", description)


def _describe_files_in_multiple_records(records, pairs):
    lines = []
    prefix = "Filename"
    for _, filename in sorted(pairs):
        lines.append("  {0} {1}".format(prefix, filename))
        prefix = " " * len(prefix)

    prefix = "Found at"
    for record in sorted(records):
        # FIXME: Show the glob that found the above files
        lines.append("  {0} {1} :: {2} :: {3}".format(prefix, *record))
        prefix = " " * len(prefix)

    return "\n".join(lines)


def _validate_prefixes(makefiles, pipeline_variant):
    logger = logging.getLogger(__name__)
    already_validated = {}
    logger.info("Validating FASTA files")
    for makefile in makefiles:
        for prefix in makefile["Genomes"].values():
            path = prefix["Path"]
            if path in already_validated:
                prefix["IndexFormat"] = already_validated[path]["IndexFormat"]
                continue

            # Must be set to a valid value, even if FASTA file does not exist
            prefix["IndexFormat"] = ".bai"

            if not os.path.exists(path):
                level = logging.ERROR if pipeline_variant == "bam" else logging.WARNING
                logger.log(level, "Reference FASTA file does not exist: %r", path)
                continue
            elif not os.path.exists(path + ".fai"):
                logger.info("Indexing FASTA at %r", path)

            try:
                contigs = FASTA.index_and_collect_contigs(path)
            except Exception as error:
                if pipeline_variant == "bam":
                    raise MakefileError(f"Error reading/indexing FASTA: {error}")
                logging.warn("Error reading/indexing FASTA: %s", error)
            else:
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
    for sample, libraries in tuple(makefile["Samples"].items()):
        for library, barcodes in tuple(libraries.items()):
            for barcode, records in tuple(barcodes.items()):
                for record in records:
                    yield sample, library, barcode, record
