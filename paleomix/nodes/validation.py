#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import os
import re

import pysam

from paleomix.node import CommandNode, Node, NodeError
from paleomix.common.fileutils import describe_files, make_dirs
from paleomix.common.utilities import chain_sorted
from paleomix.common.sequences import reverse_complement
from paleomix.tools import factory


class DetectInputDuplicationNode(Node):
    """Attempts to detect reads included multiple times as input based on the
    presence of reads with identical names AND sequences. This is compromise
    between sensitivity, specificity, and running time.

    A possible refinement would be to consider reads with the same name where
    one read is the prefix of the other (due to different amounts of trimming
    or collapsing of reads).
    """

    def __init__(self, input_files, output_file, dependencies=()):
        Node.__init__(
            self,
            description="detecting duplicate input in %s"
            % (describe_files(input_files)),
            input_files=input_files,
            output_files=output_file,
            dependencies=dependencies,
        )

    def run(self, _):
        check_bam_files(self.input_files, self._throw_node_error)

        # Everything is ok, touch the output files
        for fpath in self.output_files:
            if os.path.dirname(fpath):
                make_dirs(os.path.dirname(fpath))

            with open(fpath, "w"):
                pass

    def _throw_node_error(self, chrom, pos, records, name, seq, qual):
        message = [
            "The same read was found multiple times at position %i on %r!"
            % (pos, chrom),
            "    Name:      %r" % (name,),
            "    Sequence:  %r" % (seq,),
            "    Qualities: %r" % (qual,),
            "",
        ]

        message.append("Read was found in these BAM files:")
        for filename, records in sorted(records.items()):
            message.append("   %s in %r" % (_summarize_reads(records), filename))

        message.append("")
        message.append(
            "This indicates that the same data has been "
            "included multiple times in the project. This "
            "can be because multiple copies of the same "
            "files were used, or because one or more files "
            "contain multiple copies of the same reads. "
            "The command 'paleomix dupcheck' may be used "
            "to review the potentially duplicated data in "
            "these BAM files.\n\n"
            "If this error was a false positive, then you "
            "may execute the following command(s) to mark "
            "this test as having succeeded:"
        )

        for fpath in self.output_files:
            message.append("$ touch '%s'" % (fpath,))

        raise NodeError("\n".join(message))


class ValidateFASTQFilesNode(CommandNode):
    def __init__(
        self, input_files, output_file, offset, collapsed=False, dependencies=()
    ):
        command = factory.new(":validate_fastq")
        command.set_option("--offset", offset)
        if collapsed:
            command.set_option("--collapsed")
        command.add_multiple_values(input_files)
        command.set_kwargs(OUT_STDOUT=output_file)

        CommandNode.__init__(
            self,
            description="validating %s" % (describe_files(input_files),),
            command=command.finalize(),
            dependencies=dependencies,
        )


class ValidateFASTAFilesNode(Node):
    def __init__(self, input_files, output_file, dependencies=()):
        Node.__init__(
            self,
            description="validating %s" % (describe_files(input_files),),
            input_files=input_files,
            output_files=output_file,
            dependencies=dependencies,
        )

        assert len(self.output_files) == 1, self.output_files

    def _run(self, _config, _temp):
        for filename in self.input_files:
            check_fasta_file(filename)
        (output_file,) = self.output_files
        if os.path.dirname(output_file):
            make_dirs(os.path.dirname(output_file))
        with open(output_file, "w"):
            pass


def check_bam_files(input_files, err_func):
    handles = []
    try:
        last_pos = None
        observed_reads = collections.defaultdict(list)
        reads_iter = _open_samfiles(handles, input_files)
        references = handles[0].references

        for (record, filename) in reads_iter:
            curr_pos = (record.pos, record.tid)
            if curr_pos != last_pos:
                _process_bam_reads(observed_reads, references, last_pos, err_func)
                observed_reads.clear()
                last_pos = curr_pos

                # Stop once the trailing, unmapped reads are reached
                if record.tid == -1:
                    break

            observed_reads[record.qname].append((record, filename))
        _process_bam_reads(observed_reads, references, last_pos, err_func)
    finally:
        for handle in handles:
            handle.close()


def _open_samfiles(handles, filenames):
    sequences = []
    for filename in filenames:
        handle = pysam.AlignmentFile(filename)
        handles.append(handle)

        sequences.append(_read_samfile(handle, filename))

    return chain_sorted(*sequences, key=_key_by_tid_pos)


def _read_samfile(handle, filename):
    for record in handle:
        if record.is_unmapped and (not record.pos or record.mate_is_unmapped):
            # Ignore unmapped reads except when these are sorted
            # according to the mate position (if mapped)
            continue
        elif record.flag & 0x900:
            # Ignore supplementary / secondary alignments
            continue

        yield (record, filename)


def _process_bam_reads(observed_reads, references, position, err_func):
    for records_and_filenames in observed_reads.values():
        if len(records_and_filenames) == 1:
            # Most read-names should be obseved at most once at a position
            continue

        result = collections.defaultdict(list)
        for record, filename in records_and_filenames:
            key = (record.is_reverse, record.qname, record.seq, record.qual)
            result[key].append((filename, record))

        for (is_reverse, name, seq, qual), filenames in result.items():
            if len(filenames) == 1:
                # Two reads had same name, but different characterstics
                continue

            records = collections.defaultdict(list)
            for filename, record in filenames:
                records[filename].append(record)

            if is_reverse:
                seq = reverse_complement(seq)
                qual = qual[::-1]

            chrom = references[position[1]]
            pos = position[0]

            err_func(chrom, pos, records, name, seq, qual)


def _summarize_reads(records):
    counts = {"mate 1": 0, "mate 2": 0, "unpaired": 0}

    for record in records:
        if record.is_paired:
            if record.is_read1:
                counts["mate 1"] += 1
            elif record.is_read2:
                counts["mate 2"] += 1
            else:
                counts["unpaired"] += 1
        else:
            counts["unpaired"] += 1

    result = []
    for key, value in sorted(counts.items()):
        if value > 1:
            result.append("%i %s reads" % (value, key))
        elif value:
            result.append("%i %s read" % (value, key))

    return ", ".join(result) or "No reads"


def _key_by_tid_pos(record):
    return (record[0].tid, record[0].pos)


def check_fasta_file(filename):
    with open(filename) as handle:
        namecache = {}
        state, linelength, linelengthchanged = _NA, None, False
        for linenum, line in enumerate(handle, start=1):
            # Only \n is allowed as not all tools  handle \r
            line = line.rstrip("\n")

            if not line:
                if state in (_NA, _IN_WHITESPACE):
                    continue
                elif state == _IN_HEADER:
                    raise NodeError(
                        "Expected FASTA sequence, found empty line"
                        "\n    Filename = %r\n    Line = %r" % (filename, linenum)
                    )
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
                    raise NodeError(
                        "Empty sequences not allowed\n"
                        "    Filename = %r\n    Line = %r" % (filename, linenum - 1)
                    )
                else:
                    assert False
            else:
                if state == _NA:
                    raise NodeError(
                        "Expected FASTA header, found %r\n"
                        "    Filename = %r\n    Line = %r" % (line, filename, linenum)
                    )
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
                        raise NodeError(
                            "Lines in FASTQ files must be of same "
                            "length\n    Filename = %r\n"
                            "    Line = %r" % (filename, linenum)
                        )
                    elif linelength != len(line):
                        linelengthchanged = True
                elif state == _IN_WHITESPACE:
                    raise NodeError(
                        "Empty lines not allowed in sequences\n"
                        "    Filename = %r\n    Line = %r" % (filename, linenum)
                    )
                else:
                    assert False

        if state == _NA:
            raise NodeError(
                "File does not contain any sequences:\n"
                "    Filename = %r" % (filename,)
            )
        elif state == _IN_HEADER:
            raise NodeError(
                "File ends with an empty sequence:\n    Filename = %r" % (filename,)
            )


# Standard nucleotides + UIPAC codes
_VALID_CHARS_STR = "ACGTN" "RYSWKMBDHV"
_VALID_CHARS = frozenset(_VALID_CHARS_STR.upper() + _VALID_CHARS_STR.lower())
_NA, _IN_HEADER, _IN_SEQUENCE, _IN_WHITESPACE = range(4)


def _validate_fasta_header(filename, linenum, line, cache):
    name = line.split(" ", 1)[0][1:]
    if not name:
        raise NodeError(
            "FASTA sequence must have non-empty name\n"
            "    Filename = %r\n    Line = %r\n" % (filename, linenum)
        )
    elif not _RE_REF_NAME.match(name):
        raise NodeError(
            "Invalid name for FASTA sequence: %r\n"
            "    Filename = %r\n    Line = %r\n" % (name, filename, linenum)
        )
    elif name in cache:
        raise NodeError(
            "FASTA sequences have identical name\n"
            "    Filename = %r\n    Name = %r\n"
            "    Line 1 = %r\n    Line 2 = %r\n"
            % (filename, name, linenum, cache[name])
        )
    cache[name] = linenum


_RE_REF_NAME = re.compile("[!-()+-<>-~][!-~]*")


def _validate_fasta_line(filename, linenum, line):
    invalid_chars = frozenset(line) - _VALID_CHARS
    if invalid_chars:
        if invalid_chars == frozenset("\r"):
            raise NodeError(
                "FASTA file contains carriage-returns ('\\r')!\n"
                "Please convert file to unix format, using e.g. "
                "dos2unix.\n    Filename = %r\n" % (filename,)
            )

        raise NodeError(
            "FASTA sequence contains invalid characters\n"
            "    Filename = %r\n    Line = %r\n"
            "    Invalid characters = %r" % (filename, linenum, "".join(invalid_chars))
        )
