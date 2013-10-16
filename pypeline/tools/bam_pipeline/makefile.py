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
from pypeline.common.utilities import fill_dict
from pypeline.common.fileutils import missing_files
from pypeline.common.makefile import \
     MakefileError, \
     REQUIRED_VALUE, \
     WithoutDefaults, \
     read_makefile, \
     IsInt, \
     IsUnsignedInt, \
     IsFloat, \
     IsStr, \
     IsNone, \
     IsBoolean, \
     And, \
     Or, \
     Not, \
     ValueIn, \
     ValuesIntersect, \
     ValuesSubsetOf, \
     StringIn, \
     StringStartsWith, \
     IsListOf, \
     IsDictOf


_READ_TYPES = set(("Single", "Collapsed", "CollapsedTruncated", "Paired"))


def read_makefiles(filenames):
    makefiles = []
    for filename in filenames:
        makefile = read_makefile(filename, _VALIDATION)
        makefile = _mangle_makefile(makefile)

        makefiles.append(makefile)

    return _validate_makefiles(makefiles)


_VALID_PREFIX_NAME = \
  And(IsStr(),
      Not(Or(StringIn(["*", "Options"] + [(s + "Reads") for s in _READ_TYPES]),
             ValuesIntersect(string.whitespace, description = "contains whitespace"))))


_VALIDATION_OPTIONS = {
    # Sequencing platform, used to tag read-groups.
    "Platform" : StringIn(("CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "PACBIO"),
                          default = "ILLUMINA"),
    # Offset for quality scores in FASTQ files.
    "QualityOffset" : ValueIn((33, 64, "Solexa"),
                              default = 33),
    # Split a lane into multiple entries, one for each (pair of) file(s)
    "SplitLanesByFilenames"  : Or(IsBoolean, IsListOf(IsStr),
                                  default = True),
    # Format to use when compressing FASTQ files ("gz" or "bz2")
    "CompressionFormat" : ValueIn(("gz", "bz2"),
                                  default = "bz2"),

    "AdapterRemoval" : {
        "Version" : ValueIn(("v1.4", "v1.5+"),
                            default = "v1.5+"),
        "--pcr1"               : IsStr,
        "--pcr2"               : IsStr,
        "--maxns"              : IsUnsignedInt,
        "--minquality"         : IsUnsignedInt,
        "--trimns"             : Or(IsNone, IsBoolean),
        "--trimqualities"      : Or(IsNone, IsBoolean),
        "--collapse"           : Or(IsNone, IsBoolean),
        "--mm"                 : Or(IsFloat, IsUnsignedInt,
                                    default = 3),
        "--minlength"          : IsUnsignedInt(default = 25),
        "--minalignmentlength" : IsUnsignedInt,
        "--shift"              : IsUnsignedInt,
        "--5prime"             : IsStr,
        },

    # Which aliger/mapper to use (BWA/Bowtie2)
    "Aligners" : {
        "Program" : ValueIn(("BWA", "Bowtie2"),
                            default = "BWA"),
        "BWA" : {
            # Minimum mapping quality (PHREAD) of reads to retain
            "MinQuality" : IsUnsignedInt(default = 0),
            # Remove unmapped reads or not
            "FilterUnmappedReads" : IsBoolean(default = True),
            # Use seed region during mapping
            # Verbose name for command-line option "-l 65535"
            "UseSeed"    : IsBoolean(default = True),
            # Any number of user specific options
            StringStartsWith("-") : Or(IsListOf(IsStr, IsInt, IsFloat),
                                       Or(IsStr, IsInt, IsFloat, IsNone)),
        },
        "Bowtie2" : {
            # Minimum mapping quality (PHREAD) of reads to retain
            "MinQuality" : IsUnsignedInt(default = 0),
            # Remove unmapped reads or not
            "FilterUnmappedReads" : IsBoolean(default = True),
            # Any number of user specific options
            StringStartsWith("-") : Or(IsListOf(IsStr, IsInt, IsFloat),
                                       Or(IsStr, IsInt, IsFloat, IsNone)),
        },
    },

    # Does sample contain PCR duplicates / what to do about it.
    # True is equivalent of 'remove'.
    "PCRDuplicates"     : StringIn((True, False, 'mark', 'filter'),
                                   default = 'filter'),
    # Qualities should be rescaled using mapDamage
    "RescaleQualities"  : IsBoolean(default = False),

    # Exclude READ_TYPES from alignment/analysis
    "ExcludeReads"   : ValuesSubsetOf(_READ_TYPES,
                                      default = []),

    # Features of pipeline
    "Features"       : ValuesSubsetOf(("Raw BAM", "Realigned BAM", "Coverage", "Summary", "mapDamage", "Depths"),
                                      default = ["Realigned BAM", "Coverage", "Summary", "mapDamage", "Depths"]),
}


_VALIDATION = {
    "Options"  : _VALIDATION_OPTIONS,

    "Prefixes" : {
        _VALID_PREFIX_NAME : {
            "Path"    : IsStr(default = REQUIRED_VALUE),
            "Label"   : ValueIn(("nuclear", "mitochondrial", "chloroplast",
                                 "plasmid", "bacterial", "viral")),
            "RegionsOfInterest" : IsDictOf(IsStr, IsStr),
        },
    },

    IsStr : { # Target
        IsStr : { # Sample
            IsStr : { # Library
                IsStr     : Or(IsStr, IsDictOf(IsStr, IsStr)),

                "Options" : WithoutDefaults(_VALIDATION_OPTIONS),
            },

            "Options" : WithoutDefaults(_VALIDATION_OPTIONS),
        },

        "Options" : WithoutDefaults(_VALIDATION_OPTIONS),
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

    _split_lanes_by_filenames(makefile)

    return makefile


def _update_options(makefile):
    def _do_update_options(options, data, path):
        options = copy.deepcopy(options)
        if "Options" in data:
            if "Features" in data["Options"]:
                raise MakefileError("Features may only be specified at root level, not at %r" \
                                    % (":".join(path),))

            # Fill out missing values using those of prior levels
            options = fill_dict(destination = data.pop("Options"),
                                source      = options)

        if len(path) < 2:
            for key in data:
                if key != "Options":
                    _do_update_options(options, data[key], path + (key,))
        else:
            data["Options"] = options

    for data in makefile["Targets"].itervalues():
        _do_update_options(makefile["Options"], data, ())


def _update_prefixes(makefile):
    prefixes = {}
    for (name, values) in makefile.get("Prefixes", {}).iteritems():
        filename = values["Path"]
        if name.endswith("*"):
            records = []
            for fname in glob.glob(filename):
                name = os.path.basename(fname).split(".")[0]
                _VALID_PREFIX_NAME(("Prefixes", name), name)
                dd = copy.copy(values)
                dd["Path"] = fname

                records.append((name, dd))
            if not records:
                raise MakefileError("Did not find any matches for glob %s" % repr(filename))
        else:
            records = [(name, values)]


        for (name, record) in records:
            if name in prefixes:
                raise MakefileError("Multiple prefixes with the same name: %s" % name)

            if not record["Path"].endswith(".fasta"):
                raise MakefileError("Path for prefix %r does not end with .fasta:\n   %r" \
                                    % (name, record["Path"]))

            record["Name"]      = name
            record["Reference"] = record["Path"]
            prefixes[name]      = record

    if not prefixes:
        raise MakefileError("At least one prefix must be specified in the makefile!")
    makefile["Prefixes"] = prefixes


def _update_lanes(makefile):
    prefixes = makefile["Prefixes"]
    for samples in makefile["Targets"].itervalues():
        for libraries in samples.itervalues():
            for lanes in libraries.itervalues():
                options = lanes.pop("Options")

                for (lane, data) in lanes.iteritems():
                    lane_type = None
                    if isinstance(data, types.StringTypes):
                        lane_type = "Raw"
                    elif isinstance(data, types.DictType):
                        if all((key in _READ_TYPES) for key in data):
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
                            # Source/Current PU may differ if a lane has been split by
                            # filenames, in which case PU_src contains the original PU,
                            # while PU_cur contains a derived PU.
                            "PU_src" : barcode,
                            "PU_cur" : barcode,
                            "PG" : record["Options"]["Aligners"]["Program"],
                            "PL" : record["Options"]["Platform"].upper()}

                    record["Tags"] = tags


def _split_lanes_by_filenames(makefile):
    for (target, sample, library, barcode, record) in _iterate_over_records(makefile):
        if record["Type"] == "Raw":
            template = record["Data"]
            record["Data"] = files = paths.collect_files(template)
            split = record["Options"]["SplitLanesByFilenames"]

            if (split == True) or (isinstance(split, list) and (barcode in split)):
                if any(missing_files(file_set) for file_set in files.itervalues()):
                    raise MakefileError("Unable to split by filename for search-string '%s', did not find files" % template)
                elif any(len(v) > 1 for v in files.itervalues()):
                    template = makefile["Targets"][target][sample][library].pop(barcode)
                    keys = ("SE",) if ("SE" in files) else ("PE_1", "PE_2")

                    input_files = [files[key] for key in keys]
                    input_files_iter = itertools.izip_longest(*input_files)
                    for (index, filenames) in enumerate(input_files_iter, start = 1):
                        assert len(filenames) == len(keys)
                        assert len(filenames[0]) == len(filenames[-1])
                        new_barcode = "%s_%03i" % (barcode, index)

                        current = copy.deepcopy(template)
                        current["Data"] = dict((key, [filename]) for (key, filename) in zip(keys, filenames))
                        current["Tags"]["PU_cur"] = new_barcode

                        makefile["Targets"][target][sample][library][new_barcode] = current


def _validate_makefiles(makefiles):
    for makefile in makefiles:
        _validate_makefile_libraries(makefile)
    _validate_makefiles_duplicate_targets(makefiles)
    _validate_makefiles_duplicate_files(makefiles)

    return makefiles


def _validate_makefile_libraries(makefile):
    libraries = collections.defaultdict(set)
    for (target, sample, library, _barcode, _) in _iterate_over_records(makefile):
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
                for raw_filenames in record["Data"].itervalues():
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
    for (target, samples) in makefile["Targets"].items():
        for (sample, libraries) in samples.items():
            for (library, barcodes) in libraries.items():
                for (barcode, record) in barcodes.items():
                    yield target, sample, library, barcode, record

