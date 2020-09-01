#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Copyright (c) 2016 Mikkel Schubert <MikkelSch@gmail.com>
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
import logging
import os
import re
import tarfile

from io import TextIOWrapper

import pysam

import paleomix.yaml
from paleomix.common.formats.fasta import FASTA
from paleomix.pipelines.zonkey.common import contig_name_to_plink_name, get_sample_names

_SETTINGS_KEYS = (
    "Format",
    "Revision",
    "Plink",
    "NChroms",
    "MitoPadding",
    "SNPDistance",
)


class BAMInfo:
    def __init__(self):
        self.nuclear_contigs = {}
        self.mt_contig = None
        self.mt_length = None
        self.mt_padding = None

    @property
    def is_nuclear(self):
        return bool(self.nuclear_contigs)

    @property
    def is_mitochondrial(self):
        return bool(self.mt_contig)

    def __repr__(self):
        tmpl = "BAMInfo(nuclear=%r, mt_contig=%r, mt_length=%r, mt_padding=%r)"

        return tmpl % (
            self.nuclear_contigs,
            self.mt_contig,
            self.mt_length,
            self.mt_padding,
        )


# Format number for database file; is incremented when the format is changed.
# The 'revision' field specifies updates to the table that do not change the
# format of the database (see below).
_SUPPORTED_DB_FORMAT_MAJOR = 1
_SUPPORTED_DB_FORMAT_MINOR = 20160112

# Required columns in the 'contigs.txt' table; additional columns are ignored
_CONTIGS_TABLE_COLUMNS = frozenset(("ID", "Size", "Checksum"))
# Required columns in the 'samples.txt' table; additional non-group columns are ignored
_SAMPLES_TABLE_COLUMNS = frozenset(("ID", "Species", "Sex", "SampleID", "Publication"))
# Regular expression for parsing Group(K) columns in samples.txt
_SAMPLES_TABLE_GROUP = re.compile(r"^Group\((?P<K>.+)\)$")


class ZonkeyDBError(RuntimeError):
    pass


class ZonkeyDB:
    def __init__(self, filename):
        self.filename = filename

        log = logging.getLogger(__name__)
        log.info("Reading Zonkey database from %r" % (filename,))

        try:
            # Require that the file is not gzip / bzip2 compressed
            _check_file_compression(filename)

            with tarfile.open(filename, "r:") as tar_handle:
                log.info("Reading settings")
                self.settings = self._read_settings(tar_handle, "settings.yaml")
                log.info("Reading list of contigs")
                self.contigs = self._read_contigs_table(tar_handle, "contigs.txt")
                log.info("Reading list of samples")
                self.samples, self.groups = self._read_samples_table(
                    tar_handle, "samples.txt"
                )
                log.info("Reading mitochondrial sequences")
                self.mitochondria = self._read_mitochondria(
                    tar_handle, "mitochondria.fasta"
                )
                log.info("Reading emperical admixture distribution")
                self.simulations = self._read_simulations(tar_handle, "simulations.txt")
                log.info("Determining sample order")
                self.sample_order = self._read_sample_order(tar_handle, "genotypes.txt")
        except (OSError, tarfile.TarError) as error:
            raise ZonkeyDBError(str(error))

        self._cross_validate()

    def validate_bam(self, filename):
        """Validates a sample BAM file, checking that it is either a valid
        mitochondrial BAM (aligned against one of the referenc mt sequences),
        or that it is a valid nuclear BAM (aligned against the reference).

        Returns one of INVALID_BAMFILE, NUC_BAMFILE, and MITO_BAMFILE.
        """
        log = logging.getLogger(__name__)
        log.info("Validating BAM file %r ", filename)

        try:
            handle = pysam.AlignmentFile(filename)
        except (ValueError, IOError) as error:
            log.error("Error reading BAM: %s", error)
            return

        return self.validate_bam_handle(handle)

    def validate_bam_handle(self, handle):
        if len(get_sample_names(handle)) > 1:
            log = logging.getLogger(__name__)
            log.warning(
                "BAM read-groups specify more than one sample, "
                "but this tool treats BAMs as a single sample"
            )

        info = BAMInfo()
        if not _validate_mito_bam(self, handle, info):
            return

        if not _validate_nuclear_bam(self, handle, info):
            return

        return info

    def _cross_validate(self):
        """Cross validates tables to ensure consistency."""
        genotypes = set(self.sample_order)
        samples = set(self.samples)
        differences = (genotypes | samples) - (genotypes & samples)
        if differences:
            raise ZonkeyDBError(
                "Mismatch between samples in sample-list and "
                "genotypes table; some samples not found in "
                "both tables: %s" % (",".join(differences),)
            )

        if self.mitochondria is None:
            return

        for name, record in self.mitochondria.items():
            if name not in self.samples:
                # Ignore extra reference sequences
                meta = record.meta.upper()
                if "EXCLUDE" not in list(map(str.strip, meta.split(";"))):
                    raise ZonkeyDBError(
                        "Unexpected mitochondrial sequence: %r" % (name,)
                    )

    @classmethod
    def _read_contigs_table(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        table = cls._read_table(tar_handle, filename, _CONTIGS_TABLE_COLUMNS)
        for key, row in table.items():
            try:
                row["Size"] = int(row["Size"])
            except ValueError as error:
                raise ZonkeyDBError(
                    "Invalid size specified for sample %r in "
                    "%r: %r" % (key, filename, error)
                )

            if row["Size"] <= 0:
                raise ZonkeyDBError(
                    "Contig size must be >= 0 for %r in %r, "
                    "not %r" % (key, filename, row["Size"])
                )
        return table

    @classmethod
    def _read_samples_table(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        samples = cls._read_table(tar_handle, "samples.txt", _SAMPLES_TABLE_COLUMNS)
        if not samples:
            raise ZonkeyDBError("No samples found in genotypes table!")

        for row in samples.values():
            if row["Sex"].upper() not in ("MALE", "FEMALE", "NA"):
                raise ZonkeyDBError(
                    "Unexpected sample sex (%r); "
                    "expected 'MALE', 'FEMALE', or 'NA'" % (row["Sex"],)
                )

        group_keys = []
        for key in next(iter(samples.values())):
            match = _SAMPLES_TABLE_GROUP.match(key)
            if match is not None:
                k_value = match.groupdict()["K"]
                if not k_value.isdigit():
                    raise ZonkeyDBError(
                        "Malformed Group column name; K is " "not a number: %r" % (key,)
                    )
                elif not (2 <= int(k_value) <= 7):
                    raise ZonkeyDBError(
                        "K must be between 2 and 7, but found %r" % (key,)
                    )

                group_keys.append((key, int(k_value)))

        groups = {}
        for key, k_value in group_keys:
            group = {}
            for sample_key, sample in samples.items():
                group[sample_key] = sample.pop(key)

            group_labels = frozenset(group.values())
            if group_labels == frozenset("-"):
                continue  # Allowed for backwards compatibility
            elif "-" in group_labels:
                raise ZonkeyDBError(
                    "Not all samples column %r assignd a group" % (key,)
                )
            elif len(group_labels) != k_value:
                raise ZonkeyDBError(
                    "Expected %i groups in column %r, found %i"
                    % (k_value, key, len(group_labels))
                )

            groups[k_value] = group

        if not groups:
            raise ZonkeyDBError("No valid groups in samples.txt")

        return samples, groups

    @classmethod
    def _read_sample_order(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        handle = TextIOWrapper(tar_handle.extractfile(filename))
        header = handle.readline().rstrip("\r\n").split("\t")
        sample_order = tuple(header[-1].split(";"))

        if len(sample_order) != len(set(sample_order)):
            raise ZonkeyDBError("Duplicate sample names in %r" % (filename,))

        return sample_order

    def _read_mitochondria(self, tar_handle, filename):
        try:
            tar_handle.getmember(filename)
        except KeyError:
            # Missing MT file is allowed
            return None

        handle = TextIOWrapper(tar_handle.extractfile(filename))

        results = {}
        for record in FASTA.from_lines(handle):
            record = FASTA(
                name=record.name, meta=record.meta, sequence=record.sequence.upper()
            )

            unexpected = set(record.sequence) - set("ACGTN-")
            if unexpected:
                unexpected = ", ".join(map(repr, sorted(unexpected)))
                raise ZonkeyDBError(
                    "Unexpected nucleotide in %s; only A, C, "
                    "G, T, N, and - are allowed, not %s" % (unexpected, filename)
                )
            elif record.name in results:
                raise ZonkeyDBError(
                    "Duplicate sequence name in %s: %r" % (filename, record.name)
                )

            results[record.name] = record

        lengths = frozenset(len(record.sequence) for record in results.values())

        if not lengths:
            raise ZonkeyDBError("No mitochondrial sequences found in %r" % (filename,))
        elif len(lengths) > 2:
            lengths = tuple(sorted(lengths))
            lengths_s = "%s, and %s" % (", ".join(map(str, lengths[:-1])), lengths[-1])

            raise ZonkeyDBError(
                "At most two different sequence lengths "
                "expected for mitochondrial sequences, but "
                "found %i different lengths in %r: %s"
                % (len(lengths), filename, lengths_s)
            )
        elif len(lengths) != 1:
            # Unpadded sequences are allowed
            delta_len = max(lengths) - min(lengths)
            mito_padding = self.settings["MitoPadding"]

            if delta_len != mito_padding:
                raise ZonkeyDBError(
                    "Length difference between mitochondrial "
                    "sequences in %r does not match the "
                    "padding; expected a difference of %i bp, "
                    "but found a %i bp difference."
                    % (filename, mito_padding, delta_len)
                )

        return results

    @classmethod
    def _read_settings(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        handle = TextIOWrapper(tar_handle.extractfile(filename))

        try:
            result = paleomix.yaml.safe_load(handle)
        except paleomix.yaml.YAMLError as error:
            raise ZonkeyDBError(
                "Error reading settings file %r; %s" % (filename, error)
            )

        for key in _SETTINGS_KEYS:
            if key != "Plink":
                if not isinstance(result[key], int) or result[key] < 0:
                    raise ZonkeyDBError(
                        "Value for %r in %s must be an non-"
                        "negative integer, not %r" % (key, filename, result[key])
                    )
            elif not isinstance(result[key], str):
                raise ZonkeyDBError(
                    "Value for %r in %s must be a string, "
                    "not %r" % (key, filename, result[key])
                )

        if result["Format"] > _SUPPORTED_DB_FORMAT_MAJOR:
            raise ZonkeyDBError(
                "Database version is too old; this version of "
                "PALEOMIX supports the Zonkey DB v%i, but the "
                "database is v%i; download an updated "
                "database to continue." % (_SUPPORTED_DB_FORMAT_MAJOR, result["Format"])
            )
        elif result["Format"] < _SUPPORTED_DB_FORMAT_MAJOR:
            raise ZonkeyDBError(
                "Database version is too new; this version of "
                "PALEOMIX supports the Zonkey DB v%i, but the "
                "database is v%i; upgrade PALEOMIX to "
                "continue." % (_SUPPORTED_DB_FORMAT_MAJOR, result["Format"])
            )
        elif result["Revision"] < _SUPPORTED_DB_FORMAT_MINOR:
            raise ZonkeyDBError(
                "Database version is too old; this version of "
                "PALEOMIX supports the Zonkey DB v%i, rev. %i "
                "or newer, but the database is v%i rev. %i; "
                "please download an updated database to "
                "continue."
                % (
                    _SUPPORTED_DB_FORMAT_MAJOR,
                    _SUPPORTED_DB_FORMAT_MINOR,
                    result["Format"],
                    result["Revision"],
                )
            )

        return result

    def _read_simulations(self, tar_handle, filename):
        try:
            handle = TextIOWrapper(tar_handle.extractfile(filename))
        except KeyError:
            # Missing simulations file is allowed
            return None

        header = handle.readline().rstrip().split("\t")

        required_keys = set(
            ("NReads", "K", "Sample1", "Sample2", "HasTS", "Percentile", "Value")
        )
        missing_keys = required_keys - set(header)
        if missing_keys:
            missing_keys = ", ".join(map(repr, missing_keys))
            raise ZonkeyDBError(
                "Simulations table %r does not contain all "
                "required columns; columns %r are missing!" % (filename, missing_keys)
            )

        result = []
        for linenum, line in enumerate(handle, start=2):
            fields = line.strip().split("\t")
            if len(fields) != len(header):
                raise ZonkeyDBError(
                    "Line %i in simulations table %r, does "
                    "not contain the expected number of "
                    "columns; expected %i, but found %i!"
                    % (linenum, filename, len(header), len(fields))
                )

            row = dict(zip(header, fields))

            if row["HasTS"] not in ("TRUE", "FALSE"):
                pass

            row["HasTS"] = row["HasTS"] == "TRUE"

            for key in ("NReads", "K"):
                try:
                    row[key] = int(row[key])
                except ValueError:
                    raise ZonkeyDBError(
                        "Malformed value for column %r at "
                        "line %i in simulations table %r; "
                        "expected int, found %r" % (key, linenum, filename, row[key])
                    )

            for key in ("Percentile", "Value"):
                try:
                    row[key] = float(row[key])
                except ValueError:
                    raise ZonkeyDBError(
                        "Malformed value for column %r at "
                        "line %i in simulations table %r; "
                        "expected float, found %r" % (key, linenum, filename, row[key])
                    )

            for key in ("Sample1", "Sample2"):
                groups = frozenset(self.groups[int(row["K"])].values())

                if row[key] not in groups and row[key] != "-":
                    raise ZonkeyDBError(
                        "Invalid group in column %r in "
                        "simulations table %r: %r" % (key, filename, row[key])
                    )

            result.append(row)

        return result

    @classmethod
    def _check_required_file(cls, tar_handle, filename):
        try:
            obj = tar_handle.getmember(filename)
        except KeyError:
            raise ZonkeyDBError(
                "Database does not contain required file %r; "
                "please ensure that this is a valid Zonkey "
                "database file!" % (filename,)
            )

        if not obj.isfile():
            raise ZonkeyDBError(
                "Object %r in Zonkey database is not a "
                "file; please ensure that this is a valid "
                "Zonkey database file!" % (filename,)
            )

    @classmethod
    def _read_table(cls, tar_handle, filename, requied_columns=()):
        requied_columns = frozenset(requied_columns) | frozenset(("ID",))
        handle = TextIOWrapper(tar_handle.extractfile(filename))
        result = {}

        try:
            header = handle.readline().rstrip("\r\n").split("\t")
            if len(header) != len(set(header)):
                raise ZonkeyDBError(
                    "Table %r does contains duplicate columns!" % (filename,)
                )

            if requied_columns - set(header):
                raise ZonkeyDBError(
                    "Required columns are missign in table "
                    "%r: %s" % (filename, ", ".join())
                )

            for linenum, line in enumerate(handle):
                fields = line.rstrip("\r\n").split("\t")

                if len(fields) != len(header):
                    raise ZonkeyDBError(
                        "Error reading %r at line %i; "
                        "expected  %i columns, found %i "
                        "columns!" % (filename, linenum, len(header), len(fields))
                    )

                row = dict(zip(header, fields))
                if row["ID"] in result:
                    raise ZonkeyDBError(
                        "Duplicate IDs in %r: %s" % (filename, row["ID"])
                    )

                result[row["ID"]] = row
        finally:
            handle.close()

        return result


def _validate_mito_bam(data, handle, info):
    if data.mitochondria is None:
        # No mitochondrial data .. skip phylogeny
        return True

    references = handle.references
    min_length = min((len(record.sequence)) for record in data.mitochondria.values())
    log = logging.getLogger(__name__)

    for bam_contig, bam_length in zip(references, handle.lengths):
        if bam_contig not in data.mitochondria:
            continue

        db_sequence = data.mitochondria[bam_contig].sequence
        db_length = len(db_sequence) - db_sequence.count("-")

        if bam_length != db_length:
            log.error(
                "Length of mitochondrial contig %r (%i bp) "
                "does not match the length of the corresponding "
                "sequence in the database (%i bp)" % (bam_contig, bam_length, db_length)
            )
            return False

        filename = handle.filename.decode("utf-8")
        if not os.path.exists(filename + ".bai"):
            log.info("Indexing BAM file %r" % (filename,))
            pysam.index(filename)

        # Workaround for pysam < 0.9 returning list, >= 0.9 returning str
        for line in "".join(pysam.idxstats(filename)).split("\n"):
            line = line.strip()
            if not line:
                continue

            name, _, hits, _ = line.split("\t")
            if (name == bam_contig) and not int(hits):
                log.warning(
                    "Mitochondrial BAM (%r) does not contain "
                    "any reads aligned to contig %r; inferring an "
                    "phylogeny is not possible." % (filename, name)
                )
                return True

        info.mt_contig = bam_contig
        info.mt_length = bam_length
        info.mt_padding = len(db_sequence) - min_length

        return True
    return True


def _validate_nuclear_bam(data, handle, info):
    # Match reference panel contigs with BAM contigs; identification is done
    # by size since different repositories use different naming schemes.
    bam_contigs = collections.defaultdict(list)
    for name, length in zip(handle.references, handle.lengths):
        bam_contigs[length].append(name)

    log = logging.getLogger(__name__)
    panel_names_to_bam = {}
    for name, stats in sorted(data.contigs.items()):
        bam_contig_names = bam_contigs.get(stats["Size"], ())
        if len(bam_contig_names) == 1:
            panel_names_to_bam[name] = bam_contig_names[0]
        elif len(bam_contig_names) > 1:
            candidates = []
            for bam_name in bam_contig_names:
                if contig_name_to_plink_name(bam_name) == name:
                    candidates.append(bam_name)

            if len(candidates) == 1:
                panel_names_to_bam[name] = candidates[0]
            else:
                log.error(
                    "Multiple candidates for chr%s with size %i: %s",
                    name,
                    stats["Size"],
                    ", ".join(bam_contig_names),
                )

    if len(panel_names_to_bam) == len(data.contigs):
        info.nuclear_contigs = panel_names_to_bam
        return True
    elif panel_names_to_bam:
        log.error("Not all nuclear chromosomes found in BAM:")
        for (name, stats) in sorted(data.contigs.items()):
            is_found = "OK" if name in panel_names_to_bam else "Not found!"
            log.error("  - %s: %s" % (name, is_found))

        return False
    else:
        return True


def _check_file_compression(filename):
    with open(filename, "rb") as handle:
        header = handle.read(2)

    if header == b"\x1f\x8b":
        raise ZonkeyDBError(
            "Zonkey database is gzip compressed; please decompress to continue:\n"
            "  $ gunzip %r" % (filename,)
        )
    elif header == b"BZ":
        raise ZonkeyDBError(
            "Zonkey database is bzip2 compressed; please decompress to continue:\n"
            "  $ bunzip2 %r" % (filename,)
        )
