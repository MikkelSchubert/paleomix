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
import collections

import pypeline.tools.bam_pipeline as common


def validate_records(records):
    """Tests for possible errors in the makefile, including
     - Non-unique entries (Name ... Barcode columns)
     - Paths entered multiple times
     - Libraries with multiple samples."""
    return _validate_records_unique(records) \
        and _validate_records_libraries(records) \
        and _validate_records_paths(records)


def _validate_records_unique(records):
    paths = collections.defaultdict(list)
    for runs in records.itervalues():
        for record in runs:
            paths[common.paths.full_path(record)].append(record)

    errors = False
    for records in paths.itervalues():
        if len(records) > 1:
            errors = True
            ui.print_err("ERROR: {0} too similar records found (combination of specified fields must be unique):".format(len(records)))
            ui.print_err("\t- Name:    {Name}\n\t\t- Sample:  {Sample}\n\t\t- Library: {Library}\n\t\t- Barcode: {Barcode}\n".format(**records[0]))

    return not errors


def _validate_records_libraries(records):
    """Checks that library names are unique to each sample in a target,
    under the assumption that multiple libraries may be produced from
    a sample, but not vice versa."""
    errors = False
    for (name, runs) in sorted(records.items()):
        libraries = collections.defaultdict(set)
        for record in runs:
            libraries[record["Library"]].add(record["Sample"])

        for (library, samples) in sorted(libraries.items()):
            if len(samples) > 1:
                ui.print_err("ERROR: Multiple samples in one library:")
                ui.print_err("\t- Target:  {0}\n\t- Library: {1}\n\t- Samples: {2}\n".format(name, library, ", ".join(samples)))
                errors = True

    return not errors


def _validate_records_paths(records):
    errors = False
    for records, filenames in common.paths.collect_overlapping_files(records):
        errors = True
        descriptions = []
        for (ii, record) in enumerate(records, start = 1):
            descriptions.append("\t- Record {0}: Name: {Name},  Sample: {Sample},  Library: {Library},  Barcode: {Barcode}\n\t            Path: {Path}".format(ii, **record))
        for (ii, filename) in enumerate(filenames, start = 1):
            descriptions.append("\t- Canonical path {0}: {1}".format(ii, filename))

        ui.print_err("ERROR: Path included multiple times by one or more records:\n{0}\n".format("\n".join(descriptions)))

    return not errors


