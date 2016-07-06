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
import os
import tarfile

import pysam

import paleomix.yaml

from paleomix.common.fileutils import \
    swap_ext
from paleomix.common.formats.fasta import \
    FASTA
from paleomix.tools.zonkey.common import \
    get_sample_names, \
    contig_name_to_plink_name

from paleomix.ui import \
    print_warn, \
    print_info, \
    print_err


_SETTINGS_KEYS = ('Format', 'Revision', 'Plink', 'NChroms', 'MitoPadding',
                  'SNPDistance')


class BAMInfo(object):
    def __init__(self):
        self.nuclear = False
        self.mt_contig = None
        self.mt_length = None
        self.mt_padding = None

    @property
    def is_nuclear(self):
        return self.nuclear

    @property
    def is_mitochondrial(self):
        return bool(self.mt_contig)

    def __repr__(self):
        tmpl = "BAMInfo(nuclear=%r, mt_contig=%r, mt_length=%r, mt_padding=%r)"

        return tmpl % (self.nuclear, self.mt_contig,
                       self.mt_length, self.mt_padding)


# Format number for database file; is incremented when the format is changed.
# The 'revision' field specifies updates to the table that do not change the
# format of the database (see below).
_SUPPORTED_DB_FORMAT_MAJOR = 1
_SUPPORTED_DB_FORMAT_MINOR = 20160112

# Required columns in the 'contigs.txt' table; additional columns are ignored
_CONTIGS_TABLE_COLUMNS = frozenset(('ID', 'Size', 'Checksum'))
# Required columns in the 'samples.txt' table; additional columns are ignored
_SAMPELS_TABLE_COLUMNS = frozenset(('ID', 'Group(2)', 'Group(3)', 'Species',
                                    'Sex', 'SampleID', 'Publication'))


class ZonkeyDBError(RuntimeError):
    pass


class ZonkeyDB(object):
    def __init__(self, filename):
        self.filename = filename

        if not os.path.exists(filename):
            raise ZonkeyDBError('Database file does not exist')
        elif not tarfile.is_tarfile(filename):
            raise ZonkeyDBError('Database file is not a valid tar-file')

        print_info('Reading Zonkey database from %r ...' % (filename,))

        # Warn if file is gzip / bzip2 compressed; gives worse throughput
        _check_file_compression(filename)

        with tarfile.open(filename) as tar_handle:
            print_info('  - Reading settings ...')
            self.settings = self._read_settings(tar_handle, "settings.yaml")
            print_info('  - Reading list of contigs ...')
            self.contigs = self._read_contigs_table(tar_handle, "contigs.txt")
            print_info('  - Reading list of samples ...')
            self.samples = self._read_samples_table(tar_handle, "samples.txt")
            print_info('  - Reading mitochondrial sequences ...')
            self.mitochondria = self._read_mitochondria(tar_handle,
                                                        "mitochondria.fasta")
            print_info('  - Reading emperical admixture distribution ...')
            self.simulations = self._read_simulations(tar_handle,
                                                      "simulations.txt")
            print_info('  - Determining sample order ...')
            self.sample_order = self._read_sample_order(tar_handle,
                                                        "genotypes.txt")

        self._cross_validate()

    def validate_bam(self, filename):
        """Validates a sample BAM file, checking that it is either a valid
        mitochondrial BAM (aligned against one of the referenc mt sequences),
        or that it is a valid nuclear BAM (aligned against the reference).

        Returns one of INVALID_BAMFILE, NUC_BAMFILE, and MITO_BAMFILE.
        """
        print_info("  - Validating BAM file %r ... " % (filename,))

        try:
            handle = pysam.Samfile(filename)
        except (ValueError, IOError), error:
            print_err("Error reading BAM: %s" % (error,))
            return

        return self.validate_bam_handle(handle)

    def validate_bam_handle(self, handle):
        samples = get_sample_names(handle)
        if len(samples) > 1:
            print_warn("\nWARNING:")
            print_warn("BAM read-groups specify more than one sample, "
                       "but this tool treats BAMs as a single sample:")

            for sample in enumerate(samples, start=1):
                print_warn("    %i: %r" % sample)
            print_warn("")

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
            raise ZonkeyDBError("Mismatch between samples in sample-list and "
                                "genotypes table; some samples not found in "
                                "both tables: %s"
                                % (",".join(differences),))

        if self.mitochondria is None:
            return

        for name, record in self.mitochondria.iteritems():
            if name not in self.samples:
                # Ignore extra reference sequences
                meta = (record.meta or "").upper()
                if "EXCLUDE" not in map(str.strip, meta.split(";")):
                    raise ZonkeyDBError("Unexpected mitochondrial sequence: %r"
                                        % (name,))

    @classmethod
    def _read_contigs_table(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        table = cls._read_table(tar_handle, filename, _CONTIGS_TABLE_COLUMNS)
        for key, row in table.iteritems():
            try:
                row["Size"] = int(row["Size"])
            except ValueError, error:
                raise ZonkeyDBError("Invalid size specified for sample %r in "
                                    "%r: %r" % (key, filename, error))

            if row["Size"] <= 0:
                raise ZonkeyDBError("Contig size must be >= 0 for %r in %r, "
                                    "not %r" % (key, filename, row["Size"]))
        return table

    @classmethod
    def _read_samples_table(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        samples = cls._read_table(tar_handle, "samples.txt")
        if not samples:
            raise ZonkeyDBError("ERROR: No samples found in genotypes table!")

        for row in samples.itervalues():
            if row["Sex"].upper() not in ("MALE", "FEMALE", "NA"):
                raise ZonkeyDBError("ERROR: Unexpected sample sex (%r); "
                                    "expected 'MALE', 'FEMALE', or 'NA'"
                                    % (row["Sex"],))

        for k_groups in (2, 3):
            key = "Group(%i)" % (k_groups,)
            groups = frozenset(row[key] for row in samples.itervalues())

            if len(groups - set('-')) not in (0, k_groups):
                raise ZonkeyDBError("The %r column in the samples table must "
                                    "either contain %i ancestral groups, or "
                                    "none" % (key, k_groups))

        return samples

    @classmethod
    def _read_sample_order(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        handle = tar_handle.extractfile(filename)
        header = handle.readline().rstrip('\r\n').split('\t')
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

        handle = tar_handle.extractfile(filename)

        results = {}
        for record in FASTA.from_lines(handle):
            record = FASTA(name=record.name,
                           meta=record.meta,
                           sequence=record.sequence.upper())

            unexpected = set(record.sequence) - set("ACGTN-")
            if unexpected:
                unexpected = ", ".join(map(repr, sorted(unexpected)))
                raise ZonkeyDBError("Unexpected nucleotide in %s; only A, C, "
                                    "G, T, N, and - are allowed, not %s"
                                    % (unexpected, filename))
            elif record.name in results:
                raise ZonkeyDBError("Duplicate sequence name in %s: %r"
                                    % (filename, record.name))

            results[record.name] = record

        lengths = frozenset(len(record.sequence)
                            for record in results.itervalues())

        if not lengths:
            raise ZonkeyDBError("No mitochondrial sequences found in %r"
                                % (filename,))
        elif len(lengths) > 2:
            lengths = tuple(sorted(lengths))
            lengths_s = "%s, and %s" % (", ".join(map(str, lengths[:-1])),
                                        lengths[-1])

            raise ZonkeyDBError("At most two different sequence lengths "
                                "expected for mitochondrial sequences, but "
                                "found %i different lengths in %r: %s"
                                % (len(lengths), filename, lengths_s))
        elif len(lengths) != 1:
            # Unpadded sequences are allowed
            delta_len = max(lengths) - min(lengths)
            mito_padding = self.settings["MitoPadding"]

            if (delta_len != mito_padding):
                raise ZonkeyDBError("Length difference between mitochondrial "
                                    "sequences in %r does not match the "
                                    "padding; expected a difference of %i bp, "
                                    "but found a %i bp difference."
                                    % (filename, mito_padding, delta_len))

        return results

    @classmethod
    def _read_settings(cls, tar_handle, filename):
        cls._check_required_file(tar_handle, filename)

        handle = tar_handle.extractfile(filename)

        try:
            result = paleomix.yaml.safe_load(handle.read())
        except paleomix.yaml.YAMLError, error:
            raise ZonkeyDBError("Error reading settings file %r; %s"
                                % (filename, error))

        for key in _SETTINGS_KEYS:
            if key != "Plink":
                if not isinstance(result[key], int) or result[key] < 0:
                    raise ZonkeyDBError("Value for %r in %s must be an non-"
                                        "negative integer, not %r"
                                        % (key, filename, result[key]))
            elif not isinstance(result[key], str):
                    raise ZonkeyDBError("Value for %r in %s must be a string, "
                                        "not %r"
                                        % (key, filename, result[key]))

        if result["Format"] > _SUPPORTED_DB_FORMAT_MAJOR:
            raise ZonkeyDBError("Database version is too old; this version of "
                                "PALEOMIX supports the Zonkey DB v%i, but the "
                                "database is v%i; download an updated "
                                "database to continue."
                                % (_SUPPORTED_DB_FORMAT_MAJOR,
                                   result["Format"]))
        elif result["Format"] < _SUPPORTED_DB_FORMAT_MAJOR:
            raise ZonkeyDBError("Database version is too new; this version of "
                                "PALEOMIX supports the Zonkey DB v%i, but the "
                                "database is v%i; upgrade PALEOMIX to "
                                "continue."
                                % (_SUPPORTED_DB_FORMAT_MAJOR,
                                   result["Format"]))
        elif result["Revision"] < _SUPPORTED_DB_FORMAT_MINOR:
            raise ZonkeyDBError("Database version is too old; this version of "
                                "PALEOMIX supports the Zonkey DB v%i, rev. %i "
                                "or newer, but the database is v%i rev. %i; "
                                "please download an updated database to "
                                "continue."
                                % (_SUPPORTED_DB_FORMAT_MAJOR,
                                   _SUPPORTED_DB_FORMAT_MINOR,
                                   result["Format"],
                                   result["Revision"]))

        return result

    def _read_simulations(self, tar_handle, filename):
        try:
            handle = tar_handle.extractfile(filename)
        except KeyError:
            # Missing simulations file is allowed
            return None

        header = handle.readline().rstrip().split('\t')

        required_keys = set(('NReads', 'K', 'Sample1', 'Sample2', 'HasTS',
                             'Percentile', 'Value'))
        missing_keys = required_keys - set(header)
        if missing_keys:
            missing_keys = ', '.join(map(repr, missing_keys))
            raise ZonkeyDBError('Simulations table %r does not contain all '
                                'required columns; columns %r are missing!'
                                % (filename, missing_keys))

        result = []
        for linenum, line in enumerate(handle, start=2):
            fields = line.strip().split('\t')
            if len(fields) != len(header):
                raise ZonkeyDBError("Line %i in simulations table %r, does "
                                    "not contain the expected number of "
                                    "columns; expected %i, but found %i!"
                                    % (linenum, filename,
                                       len(header), len(fields)))

            row = dict(zip(header, fields))

            if row['HasTS'] not in ('TRUE', 'FALSE'):
                pass

            row['HasTS'] = (row['HasTS'] == 'TRUE')

            for key in ('NReads', 'K'):
                try:
                    row[key] = int(row[key])
                except ValueError:
                    raise ZonkeyDBError('Malformed value for column %r at '
                                        'line %i in simulations table %r; '
                                        'expected int, found %r'
                                        % (key, linenum, filename, row[key]))

            for key in ('Percentile', 'Value'):
                try:
                    row[key] = float(row[key])
                except ValueError:
                    raise ZonkeyDBError('Malformed value for column %r at '
                                        'line %i in simulations table %r; '
                                        'expected float, found %r'
                                        % (key, linenum, filename, row[key]))

            for key in ('Sample1', 'Sample2'):
                group_key = 'Group(%i)' % (row['K'],)
                groups = frozenset(row[group_key]
                                   for row in self.samples.itervalues())

                if row[key] not in groups and row[key] != '-':
                    raise ZonkeyDBError('Invalid group in column %r in '
                                        'simulations table %r: %r'
                                        % (key, filename, row[key]))

            result.append(row)

        return result

    @classmethod
    def _check_required_file(cls, tar_handle, filename):
        try:
            obj = tar_handle.getmember(filename)
        except KeyError:
            raise ZonkeyDBError("Database does not contain required file %r; "
                                "please ensure that this is a valid Zonkey "
                                "database file!" % (filename,))

        if not obj.isfile():
            raise ZonkeyDBError("Object %r in Zonkey database is not a "
                                "file; please ensure that this is a valid "
                                "Zonkey database file!"
                                % (filename,))

    @classmethod
    def _read_table(cls, tar_handle, filename, requied_columns=()):
        requied_columns = frozenset(requied_columns) | frozenset(("ID",))
        handle = tar_handle.extractfile(filename)
        result = {}

        try:
            header = handle.readline().rstrip('\r\n').split('\t')
            if len(header) != len(set(header)):
                raise ZonkeyDBError("Table %r does contains duplicate columns!"
                                    % (filename,))

            if requied_columns - set(header):
                raise ZonkeyDBError("Required columns are missign in table "
                                    "%r: %s" % (filename, ", ".join()))

            for linenum, line in enumerate(handle):
                fields = line.rstrip('\r\n').split('\t')

                if len(fields) != len(header):
                    raise ZonkeyDBError("Error reading %r at line %i; "
                                        "expected  %i columns, found %i "
                                        "columns!"
                                        % (filename, linenum,
                                           len(header), len(fields)))

                row = dict(zip(header, fields))
                if row["ID"] in result:
                    raise ZonkeyDBError("Duplicate IDs in %r: %s"
                                        % (filename, row["ID"]))

                result[row["ID"]] = row
        finally:
            handle.close()

        return result


def _validate_mito_bam(data, handle, info):
    if data.mitochondria is None:
        # No mitochondrial data .. skip phylogeny
        return True

    references = handle.references
    min_length = min((len(record.sequence))
                     for record in data.mitochondria.itervalues())

    for bam_contig, bam_length in zip(references, handle.lengths):
        if bam_contig not in data.mitochondria:
            continue

        db_sequence = data.mitochondria[bam_contig].sequence
        db_length = len(db_sequence) - db_sequence.count("-")

        if bam_length != db_length:
            print_err("ERROR: Length of mitochondrial contig %r (%i bp) "
                      "does not match the length of the corresponding "
                      "sequence in the database (%i bp)"
                      % (bam_contig, bam_length, db_length))
            return False

        if not os.path.exists(handle.filename + '.bai') \
                and not os.path.exists(swap_ext(handle.filename, '.bai')):
            print_info('    - Attempting to index BAM file %r!'
                       % (handle.filename,))
            pysam.index(handle.filename)

        # Workaround for pysam < 0.9 returning list, >= 0.9 returning str
        for line in "".join(pysam.idxstats(handle.filename)).split('\n'):
            line = line.strip()
            if not line:
                continue

            name, _, hits, _ = line.split('\t')
            if (name == bam_contig) and not int(hits):
                print_err("WARNING: Mitochondrial BAM (%r) does not contain "
                          "any reads aligned to contig %r; inferring an "
                          "phylogeny is not possible."
                          % (handle.filename, name))
                return True

        info.mt_contig = bam_contig
        info.mt_length = bam_length
        info.mt_padding = len(db_sequence) - min_length

        return True
    return True


def _validate_nuclear_bam(data, handle, info):
    # Check that chromosomes are of expected size; unused chroms are ignored.
    bam_contigs = dict(zip(map(contig_name_to_plink_name, handle.references),
                           handle.lengths))
    ref_contigs = data.contigs

    contigs_found = {}
    for name, stats in sorted(ref_contigs.iteritems()):
        if name not in bam_contigs:
            contigs_found[name] = False
        elif bam_contigs[name] != stats["Size"]:
            print_err("\nERROR: Chrom %r in the BAM does not match the "
                      "length specified in data file:\n"
                      "    - Expected: %i\n"
                      "    - Found: %i"
                      % (name, bam_contigs[name], stats["Size"]))

            return False
        else:
            contigs_found[name] = True

    if any(contigs_found.itervalues()):
        if not all(contigs_found.itervalues()):
            print_err("\nERROR: Not all nuclear chromosomes found in BAM:")
            for (name, stats) in sorted(ref_contigs.iteritems()):
                is_found = "Found" if contigs_found[name] else "Not found!"
                print_err("  - %s: %s" % (name, is_found))

            return False
        else:
            info.nuclear = True

    return True


def _check_file_compression(filename):
    try:
        with open(filename) as handle:
            header = handle.read(2)

            if header == "\x1f\x8b":
                print_warn('\nWARNING:\n'
                           'Zonkey database file %r is gzip compressed;\n'
                           'uncompressing the archive is recommended:\n'
                           '  $ gunzip "%s"\n' % (filename, filename))
            elif header == "BZ":
                print_warn('\nWARNING:\n'
                           'Zonkey database file %r is bzip2 compressed;\n'
                           'uncompressing the archive is recommended:\n'
                           '  $ bunzip2 "%s"\n' % (filename, filename))
    except IOError:
        # Errors are ignored at this stage
        pass
