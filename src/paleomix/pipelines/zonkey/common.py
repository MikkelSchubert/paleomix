# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import collections

from paleomix.common import yaml


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
        return chrom
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
    with open(filename) as handle:
        data = yaml.safe_load(handle)

        if not isinstance(data, dict):
            raise DBFileError("Summary file does not contain dictionary")

        results.update(data)

    return results
