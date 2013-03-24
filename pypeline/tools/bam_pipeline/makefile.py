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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
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
import copy
import glob
import types
import string
import itertools
import collections

import pypeline.tools.bam_pipeline.paths as paths
from pypeline.common.makefile import *


def read_makefiles(filenames):
    makefiles = []
    for filename in filenames:
        makefile = read_makefile(filename, _DEFAULTS, _VALIDATION)
        makefile = _mangle_makefile(makefile)

        makefiles.append(makefile)

    return _validate_makefiles(makefiles)


def _IsValidPrefixName(key, name):
    name = name.title()
    _read_types = ("Single", "Collapsed", "Paired")
    if (name in _read_types) or (name in [(s + "Reads") for s in _read_types]) or (name == "Options"):
        raise MakefileError("Prefixes cannot be named '%s', please use another name." % name)
    elif (set(name) & set(string.whitespace)):
        raise MakefileError("The label must not contain white-space, and cannot be '*': %s" % name)
    elif (set(name) & set(string.whitespace)):
        raise MakefileError("Prefix name must not contain whitespace:\n\t- Prefix: %s" % name)


def _IsValidOptions(path, value):
    validate_makefile(value, _VALIDATION["Options"], path)

    if "Features" in value:
        raise MakefileError("Features can only be specified at the top level, not at %s!" % ":".join(path))


_DEFAULTS = {
    "Options" : {
        # Sequencing platform, used to tag read-groups.
        "Platform" : "Illumina",
        # Offset for quality scores in FASTQ files.
        "QualityOffset" : 33,

        # Which aliger/mapper to use (BWA/Bowtie2)
        "Aligners" : {
            "Program" : "BWA",
            "BWA" : {
                # Minimum mapping quality (PHREAD) of reads to retain
                "MinQuality" : 0,
                # Use seed region during mapping
                # Verbose name for command-line option "-l 65535"
                "UseSeed"    : True,
            },
            "Bowtie2" : {
                # Minimum mapping quality (PHREAD) of reads to retain
                "MinQuality"  : 0,
            },
        },

        # Contains PCR duplicates, filter if true
        "PCRDuplicates"    : True,
        # Qualities should be rescaled using mapDamage
        "RescaleQualities" : False,

        # Exclude READ_TYPES from alignment/analysis
        "ExcludeReads"   : [],

        # Features of pipeline
        "Features"       : ["Realigned BAM",
                            "mapDamage",
                            "Coverage",
                            "Summary",
                            "Depths"]
    },
}


_VALIDATION = {
    "Options" : {
        # Sequencing platform, used to tag read-groups.
        "Platform" : OneOf("CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "PACBIO",  case_sensitive = False),
        # Offset for quality scores in FASTQ files.
        "QualityOffset" : OneOf(33, 64, "Solexa"),

        # Which aliger/mapper to use (BWA/Bowtie2)
        "Aligners" : {
            "Program" : OneOf("BWA", "Bowtie2"),
            "BWA" : {
                # Minimum mapping quality (PHREAD) of reads to retain
                "MinQuality" : And(IsInt, IsInRange(0, float("Inf"))),
                # Use seed region during mapping
                # Verbose name for command-line option "-l 65535"
                "UseSeed"    : IsBoolean,
            },
            "Bowtie2" : {
                # Minimum mapping quality (PHREAD) of reads to retain
                "MinQuality" : And(IsInt, IsInRange(0, float("Inf"))),
                # Any number of user specific options
                IsStrWithPrefix("-") : Or(IsListOf(IsStr, IsInt, IsFloat),
                                          Or(IsStr, IsInt, IsFloat, IsNone)),
            },
        },

        # Contains PCR duplicates, filter if true
        "PCRDuplicates"     : IsBoolean,
        # Qualities should be rescaled using mapDamage
        "RescaleQualities"  : IsBoolean,

        # Exclude READ_TYPES from alignment/analysis
        "ExcludeReads"   : AnyOf("Paired", "Single", "Collapsed"),

        # Features of pipeline
        "Features"       : AnyOf("Raw BAM", "Realigned BAM", "Coverage", "Summary", "mapDamage", "Depths"),
    },

    "Prefixes" : {
        _IsValidPrefixName : {
            "Path"    : IsStr,
            "Label"   : OneOf("nucl", "nuclear", "mito", "mitochondrial")
        },
    },

    IsStr : { # Target
        IsStr : { # Sample
            IsStr : { # Library
                IsStr     : Or(IsStr, IsDictOf(IsStr, IsStr)),
                "Options" : _IsValidOptions,
            },
        "Options" : _IsValidOptions,
        },
    },
}


def _mangle_makefile(makefile):
    makefile             = copy.deepcopy(makefile)
    makefile["Options"]  = makefile["Makefile"].pop("Options")
    makefile["Prefixes"] = makefile["Makefile"].pop("Prefixes")
    makefile["Targets"]  = makefile.pop("Makefile")

    _update_options(makefile)
    _update_prefixes(makefile)
    _update_lanes(makefile)
    _update_tags(makefile)

    return makefile


def _update_options(makefile):
    def _do_update_options(options, data, path):
        options = copy.deepcopy(options)
        if "Options" in data:
            options = apply_defaults(data.pop("Options"), options, ())

        if len(path) < 2:
            for key in data:
                if key != "Options":
                    _do_update_options(options, data[key], path + (key,))
        else:
            data["Options"] = options

    for (target, data) in makefile["Targets"].iteritems():
        _do_update_options(makefile["Options"], data, ())


def _update_prefixes(makefile):
    prefixes = {}
    for (name, values) in makefile.get("Prefixes", {}).iteritems():
        filename = values.get("Path")
        if not filename:
            raise MakefileError("Path not specified for prefix '%s'." % name)

        if name.endswith("*"):
            records = []
            for fname in glob.glob(filename):
                name = os.path.basename(fname).split(".")[0]
                _IsValidPrefixName(("Prefixes", name), name)

                records.append((name, {"Path" : fname, "Label" : values.get("Label")}))

            if not records:
                raise MakefileError("Did not find any matches for glob %s" % repr(filename))
        else:
            records = [(name, values)]


        for (name, record) in records:
            if name in prefixes:
                raise MakefileError("Multiple prefixes with the same name: %s" % name)

            reference = paths.reference_sequence(record["Path"])
            if not reference:
                raise MakefileError("""Could not find reference sequence for prefix '%s':
       Reference sequences MUST have the extensions .fasta or .fa, and
       be located at '${prefix}', '${prefix}.fa' or '${prefix}.fasta'.""" % name)

            record["Name"]      = name
            record["Reference"] = reference
            prefixes[name]      = record

    if not prefixes:
        raise MakefileError("At least one prefix must be specified in the makefile!")
    makefile["Prefixes"] = prefixes


def _update_lanes(makefile):
    prefixes = makefile["Prefixes"]
    for (target, samples) in makefile["Targets"].iteritems():
        for (sample, libraries) in samples.iteritems():
            for (library, lanes) in libraries.iteritems():
                options = lanes.pop("Options")

                for (lane, data) in lanes.iteritems():
                    lane_type = None
                    if isinstance(data, types.StringTypes):
                        lane_type = "Raw"
                    elif isinstance(data, types.DictType):
                        if all((key in ("Single", "Paired", "Collapsed")) for key in data):
                            lane_type = "Trimmed"
                        elif all((key in prefixes) for key in data):
                            lane_type = "BAMs"
                        else:
                            raise MakefileError("Error at Barcode level; keys must either be prefix-names, OR 'Paired', 'Single' or 'Collapsed'. Found: %s" \
                                                % (", ".join(data),))

                    lanes[lane] = {"Type"     : lane_type,
                                   "Data"     : data,
                                   "Options"  : options}


def _update_tags(makefile):
    for (target, samples) in makefile["Targets"].iteritems():
        for (sample, libraries) in samples.iteritems():
            for (library, barcodes) in libraries.iteritems():
                for (barcode, record) in barcodes.iteritems():
                    tags = {"Target"   : target,
                            "ID" : library,
                            "SM" : sample,
                            "LB" : library,
                            "PU" : barcode,
                            "PG" : record["Options"]["Aligners"]["Program"],
                            "PL" : record["Options"]["Platform"]}

                    record["Tags"] = tags


def _validate_makefiles(makefiles):
    for makefile in makefiles:
        _validate_makefile_libraries(makefile)
    _validate_makefiles_duplicate_targets(makefiles)
    _validate_makefiles_duplicate_files(makefiles)

    return makefiles


def _validate_makefile_libraries(makefile):
    libraries = collections.defaultdict(set)
    for (target, sample, library, barcode, record) in _iterate_over_records(makefile):
        libraries[(target, library)].add(sample)

    for ((target, library), samples) in libraries.iteritems():
        if len(samples) > 1:
            raise MakefileError("Library '%s' in target '%s' spans multiple samples: %s" \
                                    % (library, target, ", ".join(samples)))


def _validate_makefiles_duplicate_files(makefiles):
    filenames = collections.defaultdict(list)
    for makefile in makefiles:
        for (target, sample, library, barcode, record) in _iterate_over_records(makefile):
            current_filenames = []
            if record["Type"] == "Raw":
                for raw_filenames in paths.collect_files(record["Data"]).itervalues():
                    current_filenames.extend(raw_filenames)
            else:
                current_filenames.extend(record["Data"].values())

            for realpath in map(os.path.realpath, current_filenames):
                    filenames[realpath].append((target, sample, library, barcode))

    has_overlap = {}
    for (filename, records) in filenames.iteritems():
        if len(records) > 1:
            has_overlap[filename] = list(set(records))

    by_records = sorted(zip(has_overlap.values(), has_overlap.keys()))
    for (records, pairs) in itertools.groupby(by_records, lambda x: x[0]):
        descriptions = []
        for (ii, record) in enumerate(records, start = 1):
            descriptions.append("\t- Record {0}: Name: {1},  Sample: {2},  Library: {3},  Barcode: {4}".format(ii, *record))
        for (ii, (_, filename)) in enumerate(sorted(pairs), start = 1):
            descriptions.append("\t- Canonical path {0}: {1}".format(ii, filename))

        raise MakefileError("Path included multiple times by one or more records:\n{0}\n".format("\n".join(descriptions)))


def _validate_makefiles_duplicate_targets(makefiles):
    targets = set()
    for makefile in makefiles:
        for target in makefile["Targets"]:
            if target in targets:
                raise MakefileError("Target '%s' specified multiple times!" % target)
            targets.add(target)


def _iterate_over_records(makefile):
    for (target, samples) in makefile["Targets"].iteritems():
        for (sample, libraries) in samples.iteritems():
            for (library, barcodes) in libraries.iteritems():
                for (barcode, record) in barcodes.iteritems():
                    yield target, sample, library, barcode, record

