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
import itertools
import collections

from pypeline.tools.bam_pipeline.makefile import \
    DEFAULT_OPTIONS, \
    READ_TYPES, \
    MakefileError


# Valid PL values according to SAM Format Specification (v1.4-r962)
_PLATFORMS = ("CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "PACBIO")




def validate_makefiles(makefiles):
    for makefile in makefiles:
        _validate_makefile_libraries(makefile)
        _validate_makefile_options(makefile)
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


def _validate_makefile_options(makefile):
    for (target, sample, library, barcode, record) in _iterate_over_records(makefile):
        options = record["Options"]
        for (key, value) in options.iteritems():
            if key not in DEFAULT_OPTIONS:
                raise MakefileError("Unknown option: Key = '%s', value = '%s'" % (key, value))
            elif not isinstance(value, type(DEFAULT_OPTIONS[key])):
                raise MakefileError("Option value has wrong type:\n\tOption = '%s', value = '%s'\n\tExpected type = %s, found = %s" \
                                        % (key, value, DEFAULT_OPTIONS[key].__class__.__name__, value.__class__.__name__))
        
        if options["Platform"].upper() not in _PLATFORMS:
            raise MakefileError("Unknown Platform specified: '%s'" % options["Platform"])

        unknown_reads = set(options["ExcludeReads"]) - set(READ_TYPES)
        if unknown_reads:
            raise MakefileError("Unknown reads specified in option 'ExcludeReads': '%s'\n\tValid values are '%s'" \
                                    % ("', '".join(unknown_reads), "', '".join(READ_TYPES)))


def _validate_makefiles_duplicate_files(makefiles):
    filenames = collections.defaultdict(list)
    for makefile in makefiles:
        for (target, sample, library, barcode, record) in _iterate_over_records(makefile):
            for input_files in record["Reads"].get("Raw", {}).itervalues():
                for realpath in map(os.path.realpath, input_files):
                    filenames[realpath].append((target, sample, library, barcode))
            
            for filename in record["Reads"].get("Trimmed", {}).itervalues():
                realpath = os.path.realpath(filename)
                filenames[realpath].append((target, sample, library, barcode))

            for genomes in record["Reads"].get("BAM", {}).itervalues():
                for realpath in map(os.path.realpath, genomes.itervalues()):
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

