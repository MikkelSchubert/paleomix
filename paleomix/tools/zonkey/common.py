#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2016 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import collections

import paleomix.yaml

import paleomix.common.versions as versions


# Format number for database file; is incremented when the format is changed.
# The 'revision' field specifies updates to the table that do not change the
# format of the database (see below).
_SUPPORTED_DB_FORMAT = 1


RSCRIPT_VERSION = versions.Requirement(call=("Rscript", "--version"),
                                       search="version (\d+)\.(\d+)\.(\d+)",
                                       checks=versions.GE(3, 0, 0),
                                       priority=10)


class DBFileError(RuntimeError):
    pass


def get_sample_names(handle):
    samples = []
    for readgroup in handle.header.get("RG", ()):
        if "SM" in readgroup:
            samples.append(readgroup["SM"])
    return frozenset(samples)


def contig_name_to_plink_name(chrom):
    """Converts chromosome / contig name to the values expected by 'plink',
    namely a digit or X/Y, or returns None if the chromosome could not be
    identified.
    """
    if chrom.isdigit():
        return chrom.upper
    elif chrom.upper() in "XY":
        return chrom.upper()
    elif chrom.lower().startswith("chr") and chrom[3:].isdigit():
        return chrom[3:]
    elif chrom.lower() in ("chrx", "chry"):
        return chrom[3].upper()
    else:
        return None


def read_summary(filename, default="[MISSING VALUE!]"):
    results = collections.defaultdict(lambda: default)
    with open(filename) as makefile:
        string = makefile.read()
        data = paleomix.yaml.safe_load(string)

        if not isinstance(data, dict):
            raise DBFileError('Summary file does not contain dictionary')

        results.update(data)

    return results
