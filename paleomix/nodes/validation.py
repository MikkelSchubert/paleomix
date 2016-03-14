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
import collections
import io
import os
import re

import pysam

from paleomix.node import \
    Node, \
    NodeError
from paleomix.common.fileutils import \
    describe_files, \
    make_dirs
from paleomix.common.utilities import \
    chain_sorted
from paleomix.common.sequences import \
    reverse_complement

import paleomix.common.formats.fastq as fastq
import paleomix.common.procs as procs
import paleomix.common.sampling as sampling
import paleomix.tools.factory as factory


class DetectInputDuplicationNode(Node):
    """Attempts to detect reads included multiple times as input based on the
    presence of reads with identical names AND sequences. This is compromise
    between sensitivity, specificity, and running time.

    A possible refinement would be to consider reads with the same name where
    one read is the prefix of the other (due to different amounts of trimming
    or collapsing of reads).
    """

    def __init__(self, input_files, output_file, dependencies=()):
        Node.__init__(self,
                      description="<Detect Input Duplication: %s>"
                      % (describe_files(input_files)),
                      input_files=input_files,
                      output_files=output_file,
                      dependencies=dependencies)

    def run(self, _):
        handles = []
        try:
            last_pos = None
            observed_reads = collections.defaultdict(list)
            for (record, filename) in self._open_samfiles(handles, self.input_files):
                curr_pos = (record.pos, record.tid)
                if curr_pos != last_pos:
                    self._process_reads(observed_reads, self.output_files)
                    observed_reads.clear()
                    last_pos = curr_pos

                    # Stop once the trailing, unmapped reads are reached
                    if record.tid == -1:
                        break

                observed_reads[record.qname].append((record, filename))
            self._process_reads(observed_reads, self.output_files)

            # Everything is ok, touch the output files
            for fpath in self.output_files:
                make_dirs(os.path.dirname(fpath))
                with open(fpath, "w"):
                    pass
        finally:
            for handle in handles:
                handle.close()

    @classmethod
    def _open_samfiles(cls, handles, filenames):
        sequences = []
        for filename in filenames:
            handle = pysam.Samfile(filename)
            handles.append(handle)

            sequences.append(cls._read_samfile(handle, filename))

        return chain_sorted(*sequences, key=cls._key_by_tid_pos)

    @classmethod
    def _read_samfile(cls, handle, filename):
        for record in handle:
            if record.is_unmapped and (not record.pos or record.mate_is_unmapped):
                # Ignore unmapped reads except when these are sorted
                # according to the mate position (if mapped)
                continue
            elif record.flag & 0x900:
                # Ignore supplementary / secondary alignments
                continue

            yield (record, filename)

    @classmethod
    def _process_reads(cls, observed_reads, output_files):
        for records_and_filenames in observed_reads.itervalues():
            if len(records_and_filenames) == 1:
                # Most read-names should be obseved at most once at a position
                continue

            result = collections.defaultdict(list)
            for record, filename in records_and_filenames:
                key = (record.is_reverse, record.qname, record.seq, record.qual)
                result[key].append(filename)

            for (is_reverse, name, seq, qual), filenames in result.iteritems():
                if len(filenames) == 1:
                    # Two reads had same name, but different characterstics
                    continue

                filename_counts = collections.defaultdict(int)
                for filename in filenames:
                    filename_counts[filename] += 1

                if is_reverse:
                    seq = reverse_complement(seq)
                    qual = qual[::-1]

                message = ["The same read was found multiple times!",
                           "    Name:      %r" % (name,),
                           "    Sequence:  %r" % (seq,),
                           "    Qualities: %r" % (qual,),
                           ""]

                message.append("Read was found")
                for filename, count in sorted(filename_counts.iteritems()):
                    message.append("   % 2ix in %r" % (count, filename))

                message.append("")
                message.append("This indicates that the same data files have "
                               "been included multiple times in the project. "
                               "Please review the input files used in this "
                               "project, to ensure that each set of data is "
                               "included only once!\n\n"

                               "If this is not the case, then execute the "
                               "following command(s) to mark this test as "
                               "having succeeded:")

                for fpath in output_files:
                    message.append("$ touch '%s'" % (fpath,))

                raise NodeError("\n".join(message))

    @classmethod
    def _key_by_tid_pos(cls, record):
        return (record[0].tid, record[0].pos)


class ValidateFASTQFilesNode(Node):
    def __init__(self, input_files, output_file, offset, dependencies=()):
        self._offset = offset
        Node.__init__(self,
                      description="<Validate FASTQ Files: %s>"
                      % (describe_files(input_files)),
                      input_files=input_files,
                      output_files=output_file,
                      dependencies=dependencies)

    def _run(self, _config, _temp):
        check_fastq_files(self.input_files, self._offset, True)
        output_file = tuple(self.output_files)[0]
        make_dirs(os.path.dirname(output_file))
        with open(output_file, "w"):
            pass


class ValidateFASTAFilesNode(Node):
    def __init__(self, input_files, output_file, dependencies=()):
        Node.__init__(self,
                      description="<Validate FASTA Files: %s>"
                      % (describe_files(input_files)),
                      input_files=input_files,
                      output_files=output_file,
                      dependencies=dependencies)

        assert len(self.output_files) == 1, self.output_files

    def _run(self, _config, _temp):
        for filename in self.input_files:
            check_fasta_file(filename)
        output_file, = self.output_files
        make_dirs(os.path.dirname(output_file))
        with open(output_file, "w"):
            pass


def check_fastq_files(filenames, required_offset, allow_empty=False):
    for filename in filenames:
        qualities = _read_sequences(filename)
        offsets = fastq.classify_quality_strings(qualities)
        if offsets == fastq.OFFSET_BOTH:
            raise NodeError("FASTQ file contains quality scores with both "
                            "quality offsets (33 and 64); file may be "
                            "unexpected format or corrupt. Please ensure "
                            "that this file contains valid FASTQ reads from a "
                            "single source.\n    Filename = %r" % (filename,))
        elif offsets == fastq.OFFSET_MISSING:
            if allow_empty and not qualities:
                return

            raise NodeError("FASTQ file did not contain quality scores; file "
                            "may be unexpected format or corrupt. Ensure that "
                            "the file is a FASTQ file.\n    Filename = %r"
                            % (filename,))
        elif offsets not in (fastq.OFFSET_AMBIGIOUS, required_offset):
            raise NodeError("FASTQ file contains quality scores with wrong "
                            "quality score offset (%i); expected reads with "
                            "quality score offset %i. Ensure that the "
                            "'QualityOffset' specified in the makefile "
                            "corresponds to the input.\n    Filename = %s"
                            % (offsets, required_offset, filename))


def _read_sequences(filename):
    cat_call = factory.new("cat")
    cat_call.add_multiple_values((filename,))
    cat_call = cat_call.finalized_call

    cat = None
    try:
        cat = procs.open_proc(cat_call,
                              bufsize=io.DEFAULT_BUFFER_SIZE,
                              stderr=procs.PIPE,
                              stdout=procs.PIPE)
        qualities = _collect_qualities(cat.stdout, filename)

        return sampling.reservoir_sampling(qualities, 100000)
    except:
        if cat:
            cat.kill()
            cat.wait()
            cat = None
        raise
    finally:
        rc_cat = cat.wait() if cat else 0
        if rc_cat:
            message = "Error running 'paleomix cat':\n" \
                      "  Unicat return-code = %i\n\n%s" \
                      % (rc_cat, cat.stderr.read())
            raise NodeError(message)


def _collect_qualities(handle, filename):
    header = handle.readline()
    while header:
        sequence = handle.readline()
        seperator = handle.readline()
        qualities = handle.readline()

        if not header.startswith("@"):
            if header.startswith(">"):
                raise NodeError("Input file appears to be in FASTA format "
                                "(header starts with '>', expected '@'), "
                                "but only FASTQ files are supported\n"
                                "Filename = %r" % (filename,))

            raise NodeError("Input file lacks FASTQ header (expected '@', "
                            "found %r), but only FASTQ files are supported\n"
                            "    Filename = %r" % (header[:1], filename))
        elif not qualities:
            raise NodeError("Partial record found; is not 4 lines long:\n"
                            "Filename = %r\n    Record = '%s'"
                            % (filename, header.rstrip()))
        elif not seperator.startswith("+"):
            raise NodeError("Input file lacks FASTQ seperator (expected '+', "
                            "found %r), but only FASTQ files are supported\n"
                            "    Filename = %r" % (seperator[:1], filename))
        elif len(sequence) != len(qualities):
            raise NodeError("Input file contains malformed FASTQ records; "
                            "length of sequence / qualities are not the "
                            "same.\n    Filename = %r\n    Record = '%s'"
                            % (filename, header.rstrip()))

        yield qualities
        header = handle.readline()


def check_fasta_file(filename):
    with open(filename) as handle:
        namecache = {}
        state, linelength, linelengthchanged = _NA, None, False
        for linenum, line in enumerate(handle, start=1):
            # Only \n is allowed as not all tools (e.g. GATK) handle \r
            line = line.rstrip('\n')

            if not line:
                if state in (_NA, _IN_WHITESPACE):
                    continue
                elif state == _IN_HEADER:
                    raise NodeError("Expected FASTA sequence, found empty line"
                                    "\n    Filename = %r\n    Line = %r"
                                    % (filename, linenum))
                elif state == _IN_SEQUENCE:
                    state = _IN_WHITESPACE
                else:
                    assert False
            elif line.startswith(">"):
                if state in (_NA, _IN_SEQUENCE, _IN_WHITESPACE):
                    _validate_fasta_header(filename, linenum, line, namecache)
                    state = _IN_HEADER
                    linelength = None
                    linelengthchanged = False
                elif state == _IN_HEADER:
                    raise NodeError("Empty sequences not allowed\n"
                                    "    Filename = %r\n    Line = %r"
                                    % (filename, linenum - 1))
                else:
                    assert False
            else:
                if state == _NA:
                    raise NodeError("Expected FASTA header, found %r\n"
                                    "    Filename = %r\n    Line = %r"
                                    % (line, filename, linenum))
                elif state == _IN_HEADER:
                    _validate_fasta_line(filename, linenum, line)
                    linelength = len(line)
                    state = _IN_SEQUENCE
                elif state == _IN_SEQUENCE:
                    _validate_fasta_line(filename, linenum, line)
                    # If the length has changed, then that line must be the
                    # last line in the record, which may be shorter due to the
                    # sequence length. This is because the FAI index format
                    # expects that each line has the same length.
                    if linelengthchanged or (linelength < len(line)):
                        raise NodeError("Lines in FASTQ files must be of same "
                                        "length\n    Filename = %r\n"
                                        "    Line = %r" % (filename, linenum))
                    elif linelength != len(line):
                        linelengthchanged = True
                elif state == _IN_WHITESPACE:
                    raise NodeError("Empty lines not allowed in sequences\n"
                                    "    Filename = %r\n    Line = %r"
                                    % (filename, linenum))
                else:
                    assert False

        if state in (_NA, _IN_HEADER):
            raise NodeError("File does not contain any sequences"
                            "    Filename = %r" % (filename, ))

# Standard nucleotides + UIPAC codes
_VALID_CHARS_STR = "ACGTN" "RYSWKMBDHV"
_VALID_CHARS = frozenset(_VALID_CHARS_STR.upper() + _VALID_CHARS_STR.lower())
_NA, _IN_HEADER, _IN_SEQUENCE, _IN_WHITESPACE = range(4)


def _validate_fasta_header(filename, linenum, line, cache):
    name = line.split(" ", 1)[0][1:]
    if not name:
        raise NodeError("FASTA sequence must have non-empty name\n"
                        "    Filename = %r\n    Line = %r\n"
                        % (filename, linenum))
    elif not _RE_REF_NAME.match(name):
        raise NodeError("Invalid name for FASTA sequence: %r\n"
                        "    Filename = %r\n    Line = %r\n"
                        % (name, filename, linenum))
    elif name in cache:
        raise NodeError("FASTA sequences have identical name\n"
                        "    Filename = %r\n    Name = %r\n"
                        "    Line 1 = %r\n    Line 2 = %r\n"
                        % (filename, name, linenum, cache[name]))
    cache[name] = linenum
_RE_REF_NAME = re.compile("[!-()+-<>-~][!-~]*")


def _validate_fasta_line(filename, linenum, line):
    invalid_chars = frozenset(line) - _VALID_CHARS
    if invalid_chars:
        if invalid_chars == frozenset('\r'):
            raise NodeError("FASTA file contains carriage-returns ('\\r')!\n"
                            "Please convert file to unix format, using e.g. "
                            "dos2unix.\n    Filename = %r\n" % (filename,))

        raise NodeError("FASTA sequence contains invalid characters\n"
                        "    Filename = %r\n    Line = %r\n"
                        "    Invalid characters = %r"
                        % (filename, linenum, "".join(invalid_chars)))
