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

from collections import defaultdict
from typing import (
    IO,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    cast,
)

from paleomix.common.fileutils import PathTypes, open_rt
from paleomix.common.formats.fasta import FASTA, FASTAError
from paleomix.common.sequences import NT_CODES, encode_genotype, split
from paleomix.common.utilities import safe_coerce_to_frozenset


class MSAError(FASTAError):
    pass


class MSA(FrozenSet[FASTA]):
    """Represents a Multiple Sequence Alignment of FASTA records."""

    def __new__(cls, sequences: Iterable[FASTA]) -> "MSA":
        records: List[FASTA] = []
        names: Set[str] = set()
        for record in sequences:
            if record.name in names:
                raise MSAError(
                    "Duplicate name found in FASTA records: %r" % record.name
                )
            records.append(record)
            names.add(record.name)

        if not records:
            raise MSAError("MSA does not contain any sequences")

        instance = super(MSA, cls).__new__(cls, records)  # type: ignore
        instance = cast(MSA, instance)

        MSA.validate(instance)
        return instance

    def seqlen(self) -> int:
        """Retrurns the length of the sequences in the MSA."""
        return len(next(iter(self)).sequence)

    def exclude(self, names: Iterable[str]) -> "MSA":
        """Builds a new MSA that excludes the named set of records."""
        _, excluded, _ = self._group(names)
        return MSA(excluded)

    def select(self, names: Iterable[str]) -> "MSA":
        """Builds a new MSA that includes only the named set of records."""
        included, _, _ = self._group(names)
        return MSA(included)

    def reduce(self) -> Optional["MSA"]:
        columns: List[Sequence[str]] = []
        uncalled = frozenset("Nn-")
        for column in zip(*(record.sequence for record in self)):
            if frozenset(column) - uncalled:
                columns.append(column)

        if not columns:
            return None

        records: List[FASTA] = []
        for record, sequence in zip(self, zip(*columns)):
            records.append(FASTA(record.name, record.meta, "".join(sequence)))

        return MSA(records)

    def filter_singletons(self, to_filter: str, filter_using: Iterable[str]) -> "MSA":
        included, excluded_, to_filter_ = self._group(filter_using, to_filter)
        if to_filter_ is None:
            raise KeyError(to_filter)

        sequence = list(to_filter_.sequence)
        sequences = [record.sequence.upper() for record in included]
        for index, nts in enumerate(zip(*sequences)):
            current_nt = sequence[index].upper()
            if current_nt in "N-":
                continue

            allowed_nts = set()  # type: Set[str]
            for allowed_nt in nts:
                if allowed_nt not in "N-":
                    allowed_nts.update(NT_CODES[allowed_nt])
            filtered_nts = frozenset(NT_CODES[current_nt]) & allowed_nts

            if not filtered_nts:
                filtered_nts = frozenset("N")

            genotype = encode_genotype(filtered_nts)
            if genotype != current_nt:
                sequence[index] = genotype.lower()
        new_record = FASTA(to_filter_.name, to_filter_.meta, "".join(sequence))

        return MSA([new_record] + included + excluded_)

    def split(self, split_by: str = "123") -> Dict[str, "MSA"]:
        """Splits a MSA and returns a dictionary of keys to MSAs,
        using the keys in the 'split_by' parameter at the top
        level. See also paleomix.common.sequences.split."""
        self.validate(self)
        if not split_by:
            raise TypeError("No partitions to split by specified")

        results = {key: set() for key in split_by}  # type: Dict[str, Set[FASTA]]
        for record in self:
            for key, partition in split(record.sequence, split_by).items():
                results[key].add(FASTA(record.name, None, partition))

        return {key: MSA(value) for key, value in results.items()}

    @classmethod
    def join(cls, *msas: "MSA") -> "MSA":
        """Merge multiple MSAs into a single MSA, by concatenating sequences in
        the order of the passed MSAs. Sequences are joined by name, and all MSAs
        must therefore contain the same set of sequence names. Meta information
        is not preserved."""
        cls.validate(*msas)

        merged = defaultdict(list)  # type: Dict[str, List[str]]
        for msa in msas:
            for record in msa:
                merged[record.name].append(record.sequence)

        sequences = []  # type: List[FASTA]
        for name, sequence in merged.items():
            sequences.append(FASTA(name, None, "".join(sequence)))
        return MSA(sequences)

    @classmethod
    def from_lines(cls, lines: Iterable[str]) -> "MSA":
        """Parses a MSA from a file/list of lines, and returns a dictionary
        of names to sequences. If read_meta is True, meta information included
        after the first space in header of each sequence:
          >NAME META-INFORMATION
          SEQUENCE
        As suggested above, sequences are expected to be in FASTA format."""
        return MSA(FASTA.from_lines(lines))

    @classmethod
    def from_file(cls, filename: PathTypes) -> "MSA":
        """Reads a MSA from the specified filename. The file may
        be uncompressed, gzipped or bzipped. See also 'MSA.from_lines'."""
        with open_rt(filename) as handle:
            return MSA.from_lines(handle)

    def to_file(self, fileobj: IO[str]) -> None:
        for fst in sorted(self):
            fst.write(fileobj)

    @classmethod
    def validate(cls, *msas: "MSA") -> None:
        """Validates one or more MSAs, requiring:
        1. Identical sets of sequence names across all MSAs.
        2. That all names are non-empty strings.
        3. That all sequences are of the same length (per MSA).
        4. That no empty MSA (no sequences) are specified."""
        if not msas:
            raise TypeError("No MSAs given as arguments")

        seqs_all = msas[0].names()
        seqs_common = set(seqs_all)
        for msa in msas:
            if len(set(len(record.sequence) for record in msa)) != 1:
                raise MSAError("MSA contains sequences of differing lengths")

            seqs_all.update(msa.names())
            seqs_common &= set(msa.names())

        if seqs_all != seqs_common:
            raise MSAError(
                "Some sequences not found in all MSAs: '%s'"
                % ("', '".join(seqs_all - seqs_common),)
            )

    def __repr__(self) -> str:
        def _fasta_to_str(fst: FASTA) -> str:
            return "FASTA(%r, %r, %r)" % (fst.name, fst.meta, fst.sequence)

        return "MSA(%s)" % (", ".join(map(_fasta_to_str, sorted(self))))

    def names(self) -> Set[str]:
        return set(record.name for record in self)

    def _group(
        self,
        selection: Iterable[str],
        extra: Optional[str] = None,
    ) -> Tuple[List[FASTA], List[FASTA], Optional[FASTA]]:
        selection = safe_coerce_to_frozenset(selection)
        if extra in selection:
            raise MSAError("Key used for multiple selections: %r" % extra)
        elif not selection:
            raise ValueError("No FASTA names given")

        missing_keys = selection - self.names()
        if missing_keys:
            raise KeyError("Key(s) not found: %r" % (", ".join(map(str, missing_keys))))

        included = []  # type: List[FASTA]
        excluded = []  # type: List[FASTA]
        other: Optional[FASTA] = None
        for record in self:
            if record.name in selection:
                included.append(record)
            elif record.name != extra:
                excluded.append(record)
            else:
                other = record

        return included, excluded, other
