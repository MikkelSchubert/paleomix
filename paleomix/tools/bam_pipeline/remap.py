#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
from __future__ import print_function

import os
import sys
import bz2
import pysam

from optparse import OptionParser

import paleomix.tools.bam_pipeline.config as bam_cfg
import paleomix.tools.bam_pipeline.mkfile as bam_mkfile
import paleomix.tools.bam_pipeline.pipeline as bam_pipeline
from paleomix.common.fileutils import make_dirs, reroot_path, add_postfix
from paleomix.common.sequences import reverse_complement
from paleomix.common.utilities import set_in, get_in


# Indentation used when printing makefile
_INDENTATION = " " * 4


def samefile(fname_a, fname_b):
    if not (os.path.exists(fname_a) and os.path.exists(fname_b)):
        return False
    return os.path.samefile(fname_a, fname_b)


class ReadSink(object):
    # Cache of opened file-handles; used to ensure that lanes containing both
    # PE (potentialy with orphaned reads) and SE reads are all collected.
    _cache = {}

    def __init__(self, filename, handle):
        """See ReadSink.open"""
        self.filename = filename
        self._handle = handle

    def close(self):
        self._handle.close()

    def write_records(self, records):
        for record in records:
            seq = record.seq
            qual = record.qual
            if record.is_reverse:
                seq = reverse_complement(seq)
                qual = qual[::-1]

            assert len(qual) == len(seq), record.qname
            self._handle.write("@%s\n" % (record.qname,))
            self._handle.write("%s\n" % (seq,))
            self._handle.write("+\n")
            self._handle.write("%s\n" % (qual,))

    @classmethod
    def get_filename(cls, destination, prefix):
        return os.path.join(destination, "%s.fastq.bz2" % (prefix,))

    @classmethod
    def open(cls, prefix, filename):
        if filename not in cls._cache:
            handle = bz2.BZ2File(os.path.join(prefix, filename), "w")
            cls._cache[filename] = ReadSink(filename, handle)

        return cls._cache[filename]

    @classmethod
    def close_all_sinks(cls):
        for handle in cls._cache.itervalues():
            handle.close()
        cls._cache.clear()


class PEReadSink(ReadSink):
    def __init__(self, prefix, destination):
        ReadSink.__init__(self, self.get_filename(destination, "paired.{Pair}"), None)
        self._sink_se = ReadSink.open(prefix, self.get_filename(destination, "singleton"))
        self._sink_pe_1 = ReadSink.open(prefix, self.get_filename(destination, "paired.1"))
        self._sink_pe_2 = ReadSink.open(prefix, self.get_filename(destination, "paired.2"))

    def write_records(self, records):
        record_cache = {}
        for record in records:
            num = 0
            if record.is_read1:
                num = 1
            elif record.is_read2:
                num = 2
            set_in(record_cache, (record.qname, num), record)

        for pair in record_cache.itervalues():
            # Only write complete pairs
            if (1 in pair) and (2 in pair):
                self._sink_pe_1.write_records([pair.pop(1)])
                self._sink_pe_2.write_records([pair.pop(2)])

            # Any orphan files are written to the SE sink
            for record in pair.itervalues():
                self._sink_se.write_records([record])

    @classmethod
    def open(cls, prefix, destination):
        return PEReadSink(prefix, destination)


def convert_reads(config, destination, record, sink_cache):
    # Source name is used, to re-merge split lanes
    name = record.tags.get("PU_src")
    destination = os.path.join(destination, name)
    make_dirs(os.path.join(config.destination, destination))

    def _open_se_sink(reads_type):
        key = (name, reads_type)
        if not get_in(sink_cache, key):
            filename = ReadSink.get_filename(destination, reads_type.lower())
            set_in(sink_cache, key, ReadSink.open(config.destination, filename))
        return key

    for (reads_type, bam_files) in record.bams.iteritems():
        # Processed reads are pre-aligned BAMs which have been cleaned up
        if reads_type in ("Paired", "Processed"):
            # Record "Single" reads; these may result from orphan SE reads
            _open_se_sink("Singleton")

            key = (name, "Paired")
            if not get_in(sink_cache, key):
                set_in(sink_cache, key, PEReadSink.open(config.destination,
                                                        destination))
        else:
            key = _open_se_sink(reads_type)

        sink = get_in(sink_cache, key)
        for filename in bam_files:
            print("%sProcessing file %r" % (_INDENTATION * 4, filename))
            with pysam.Samfile(filename) as handle:
                def _keep_record(record):
                    return (record.qual >= config.min_quality) and \
                        (len(record.seq) >= config.min_length)

                sink.write_records(record for record in handle
                                   if _keep_record(record))


def parse_options(argv):
    tmpl = "%s <Prefix> <Makefile> [<Makefile>, ...]"

    parser = OptionParser(tmpl % "bam_pipeline remap")
    parser.add_option("--destination",
                      default="remapping", dest="destination",
                      help="Destination for resulting files [%default]")
    parser.add_option("--output-name-postfix",
                      default="_remapping", dest="postfix",
                      help="Postfix added to filenames/target names of "
                      "generated files [%default]")
    parser.add_option("--min-quality", default=0, type=int,
                      help="Minimum quality of hits to include in output "
                           "[%default]")
    parser.add_option("--min-length", default=0, type=int,
                      help="Minimum length of hits to include in output "
                           "[%default]")

    config, args = parser.parse_args(argv)
    if (len(args) < 2):
        parser.print_usage()
        return None, None

    config.prefix = args[0]

    return config, args[1:]


def main(argv):
    config, args = parse_options(argv)
    if config is None:
        return 1

    # Get default options for bam_pipeline
    bam_config, _ = bam_cfg.parse_config(args, "bam")
    makefiles = bam_pipeline.read_makefiles(bam_config, args)
    # Build .fai files for reference .fasta files
    bam_pipeline.index_references(bam_config, makefiles)

    for makefile in makefiles:
        mkfile_fname = makefile["Statistics"]["Filename"]
        bam_config.destination = os.path.dirname(mkfile_fname)
        tasks = bam_pipeline.build_pipeline_full(bam_config, makefile,
                                                 return_nodes=False)

        make_dirs(config.destination)
        makefile_name = add_postfix(makefile["Statistics"]["Filename"],
                                    config.postfix)
        makefile_path = reroot_path(config.destination, makefile_name)
        if samefile(makefile["Statistics"]["Filename"], makefile_path):
            sys.stderr.write("ERROR: Would overwrite source makefile at %r\n" % (makefile_path,))
            sys.stderr.write("       Please set --destination and/or --output-name-postfix\n")
            sys.stderr.write("       before continuing.\n")
            return 1

        print("Writing makefile", makefile_path)

        found_prefix = False
        for prefix in makefile["Prefixes"]:
            if prefix != config.prefix:
                print("%sSkipping %s" % (_INDENTATION, prefix))
            else:
                found_prefix = True

        if not found_prefix:
            sys.stderr.write("\nERROR:\n")
            sys.stderr.write("Could not find prefix %r in %r! Aborting ...\n"
                             % (config.prefix, mkfile_fname))
            return 1

        with open(makefile_path, "w") as makefile_handle:
            template = bam_mkfile.build_makefile(add_sample_tmpl=False)
            makefile_handle.write(template)
            makefile_handle.write("\n" * 3)

            for target in tasks:
                target_name = add_postfix(target.name, config.postfix)
                print("%sTarget: %s -> %s" % (_INDENTATION,
                                              target.name,
                                              target_name))

                makefile_handle.write('%s"%s":\n' % (_INDENTATION * 0,
                                                     target_name))
                for prefix in target.prefixes:
                    if prefix.name != config.prefix:
                        continue

                    for sample in prefix.samples:
                        print("%sSample: %s" % (_INDENTATION * 2, sample.name))

                        makefile_handle.write('%s"%s":\n' % (_INDENTATION * 1,
                                                             sample.name))

                        for library in sample.libraries:
                            print("%sLibrary: %s" % (_INDENTATION * 3,
                                                     library.name))
                            makefile_handle.write('%s"%s":\n'
                                                  % (_INDENTATION * 2,
                                                     library.name))

                            sink_cache = {}
                            destination = os.path.join(target_name,
                                                       "reads",
                                                       sample.name,
                                                       library.name)

                            for lane in library.lanes:
                                convert_reads(config, destination, lane, sink_cache)
                            ReadSink.close_all_sinks()

                            for lane_name in sorted(sink_cache):
                                makefile_handle.write('%s"%s":\n' % (_INDENTATION * 3, lane_name))
                                for (reads_type, sink) in sorted(sink_cache[lane_name].items()):
                                    makefile_handle.write('%s%s "%s"\n'
                                                          % (_INDENTATION * 4,
                                                             ("%s:" % (reads_type,)).ljust(20),
                                                             sink.filename))
                                makefile_handle.write("\n")
        print("\tDone ...")
        print()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
