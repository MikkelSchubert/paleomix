# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import os
import sys
from collections.abc import Iterable, Iterator
from typing import IO

import pysam

from paleomix.common.fileutils import fspath, open_rt
from paleomix.common.formats._common import FormatError
from paleomix.common.utilities import Immutable, TotallyOrdered, fragment, split_before


class FASTAError(FormatError):
    pass


class FASTA(TotallyOrdered, Immutable):
    name: str
    meta: str | None
    sequence: str

    def __init__(self, name: str, meta: str | None, sequence: str) -> None:
        if not (name and isinstance(name, str)):
            raise FASTAError("FASTA name must be a non-empty string")
        elif not (isinstance(meta, str) or (meta is None)):
            raise FASTAError("FASTA meta must be a string, or None")
        elif not isinstance(sequence, str):
            raise FASTAError("FASTA sequence must be a string")

        Immutable.__init__(self, name=name, meta=meta or "", sequence=sequence)

    def write(self, fileobj: IO[str] = sys.stdout) -> None:
        """Writes a FASTA record to fileobj, wrapping sequences at 60 chars"""
        name = self.name
        if self.meta:
            name = f"{name} {self.meta}"

        fileobj.write(">{}\n{}\n".format(name, "\n".join(fragment(60, self.sequence))))

    @classmethod
    def from_lines(cls, lines: Iterable[str]) -> Iterator[FASTA]:
        """Parses FASTA sequences found in a sequence of lines, and returns
        a tuple for each FASTA record: ((name, meta-information), sequence)
        No assumptions are made about the line-lengths."""
        lines = (line.rstrip() for line in lines)
        for record in split_before(lines, lambda v: v.startswith(">")):
            name = record[0]
            if (not name.startswith(">")) or (len(name) == 1):
                raise FASTAError("Unnamed FASTA record")
            elif len(record) == 1:
                raise FASTAError(f"FASTA record does not contain sequence: {name[1:]}")

            # Split out any meta information
            name_and_meta = name[1:].split(None, 1)
            if len(name_and_meta) < 2:
                name_and_meta.append("")
            name, meta = name_and_meta

            yield FASTA(name=name, meta=meta, sequence="".join(record[1:]))

    @classmethod
    def from_file(cls, filename: str) -> Iterator[FASTA]:
        """Reads an unindexed FASTA file, returning a sequence of
        tuples containing the name and sequence of each entry in
        the file. The FASTA file may be GZIP/BZ2 compressed."""
        with open_rt(filename) as fasta_file:
            yield from FASTA.from_lines(fasta_file)

    @classmethod
    def index_and_collect_contigs(cls, filename: os.PathLike[str]) -> dict[str, int]:
        """Creates an index (.fai; if it does not already exist) for a FASTA file using
        pysam, and returns a dictionary of {contig: length} listed in that file.
        """
        # If an index does not already exist, then it is created automatically
        with pysam.FastaFile(fspath(filename)) as handle:
            return dict(zip(handle.references, handle.lengths))

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, FASTA):
            return NotImplemented

        return (self.name, self.meta, self.sequence) < (
            other.name,
            other.meta,
            other.sequence,
        )

    def __hash__(self) -> int:
        return hash((self.name, self.meta, self.sequence))

    def __repr__(self) -> str:
        return f"FASTA({self.name!r}, {self.meta!r}, {self.sequence!r})"
