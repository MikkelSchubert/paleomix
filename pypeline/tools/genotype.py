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
"""Wrapper around "samtools mpileup | bcftools view", to improve performance
when genotyping sparse regions (e.g. sets of exons / RNAs / or similar), and
allow transparent multi-threading. Alternatively, only "samtools mpileup" may
be called, in order to generate pileups for a set of regions.

There are 3 main motivations for this wrapper:
    1. Current versions of SAMTools read the full contents of the input BAM,
       when a BED file of regions is specified, even if these regions cover
       only a fraction of sites. This can be somewhat mitigated by ALSO
       specifying a region using -r, which fetches just that region, but this
       does not scale well for thousands of individual regions.
    2. It provides transparent parallelization, allowing any set of bed regions
       to be split and processed in parallel.
"""
import os
import sys
import shutil
import signal
import argparse
import traceback
import subprocess
import multiprocessing

import pysam

from pypeline.common.bedtools import \
    read_bed_file, \
    sort_bed_by_bamfile

from pypeline.nodes.samtools import \
    samtools_compatible_wbu_mode


class BatchError(RuntimeError):
    pass


# Size of smallest block in (linear) BAM index (= 2 << 14)
_BAM_BLOCK_SIZE = 16384


###############################################################################
###############################################################################

def popen(call, *args, **kwargs):
    """Equivalent to subprocess.Popen, but records the system call as a tuple
    assigned to the .call property of the Popen object.
    """
    proc = subprocess.Popen(call, *args, **kwargs)
    proc.call = tuple(call)
    return proc


###############################################################################
###############################################################################
# CLI functions

def build_call(call, args, positional, new_args):
    call = list(call)
    args = dict(args)

    for new_arg in new_args:
        key, value = new_arg, None
        if "=" in new_arg:
            key, value = new_arg.split("=", 1)
        args[key] = value

    for (key, value) in sorted(args.iteritems()):
        call.append(key)
        if value is not None:
            call.append(value)

    call.extend(positional)

    return call


###############################################################################
###############################################################################
# Common functions

def cleanup_batch(setup):
    print "Cleaning up batch ..."
    for handle in setup["handles"].itervalues():
        handle.close()

    for proc in setup["procs"].itervalues():
        if proc.poll() is None:
            proc.terminate()
            proc.wait()

    for filename in setup["temp_files"].itervalues():
        print "Removing temporary file %r" % (filename,)
        os.remove(filename)


def write_bed_file(prefix, regions):
    fpath = prefix + ".bed"
    with open(fpath, "w") as bed_handle:
        for (contig, start, end) in regions:
            bed_handle.write("%s\t%i\t%i\n" % (contig, start, end))
        bed_handle.flush()
    return fpath


def make_bam_pipe(prefix):
    fpath_pipe = prefix + ".pipe"
    if os.path.exists(fpath_pipe):
        os.remove(fpath_pipe)
    os.mkfifo(fpath_pipe)
    return fpath_pipe


def setup_basic_batch(args, regions, prefix, func):
    setup = {"files": {},
             "temp_files": {},
             "procs": {},
             "handles": {}}

    try:
        setup["files"]["bed"] = write_bed_file(prefix, regions)
        setup["temp_files"]["bed"] = setup["files"]["bed"]

        setup["files"]["pipe"] = make_bam_pipe(prefix)
        setup["temp_files"]["pipe"] = setup["files"]["pipe"]

        setup["handles"]["outfile"] = open(prefix, "w")
        zip_proc = popen(["bgzip"],
                         stdin=func(setup),
                         stdout=setup["handles"]["outfile"],
                         close_fds=True)

        setup["procs"]["gzip"] = zip_proc

        write_mode = samtools_compatible_wbu_mode()
        setup["handles"]["bam_in"] = pysam.Samfile(args.bamfile)
        setup["handles"]["bam_out"] = \
            pysam.Samfile(setup["files"]["pipe"], write_mode,
                          template=setup["handles"]["bam_in"])

        return setup
    except:
        traceback.print_exc()
        cleanup_batch(setup)
        raise


###############################################################################
###############################################################################
# Pileup batch generation

def setup_mpileup_batch(args, regions, prefix):
    def _create_mpileup_proc(setup):
        mpileup_args = {"-l": setup["files"]["bed"]}
        call = build_call(call=("samtools", "mpileup"),
                          args=mpileup_args,
                          new_args=args.mpileup_argument,
                          positional=(setup["files"]["pipe"],))

        print "Running 'samtools mpileup': %s" % (" ".join(call))
        procs = setup["procs"]
        procs["mpileup"] = popen(call,
                                 stdout=subprocess.PIPE,
                                 close_fds=True)

        return procs["mpileup"].stdout

    return setup_basic_batch(args, regions, prefix, _create_mpileup_proc)


###############################################################################
###############################################################################
# Genotyping batch generation

def setup_genotyping_batch(args, regions, prefix):
    def _create_genotyping_proc(setup):
        mpileup_args = {"-u": None,
                        "-l": setup["files"]["bed"]}
        mpileup_call = build_call(call=("samtools", "mpileup"),
                                  args=mpileup_args,
                                  new_args=args.mpileup_argument,
                                  positional=(setup["files"]["pipe"],))

        print "Running 'samtools mpileup': %s" % (" ".join(mpileup_call))
        procs = setup["procs"]
        procs["mpileup"] = popen(mpileup_call,
                                 stdout=subprocess.PIPE,
                                 close_fds=True)

        bcftools_call = build_call(call=("bcftools", "view"),
                                   args={},
                                   new_args=args.bcftools_argument,
                                   positional=("-",))

        print "Running 'bcftools call': %s" % (" ".join(bcftools_call))
        procs["bcftools"] = popen(bcftools_call,
                                  stdin=procs["mpileup"].stdout,
                                  stdout=subprocess.PIPE,
                                  close_fds=True)

        return procs["bcftools"].stdout

    return setup_basic_batch(args, regions, prefix, _create_genotyping_proc)


###############################################################################
###############################################################################

def setup_batch(args, regions, filename):
    """Setup a batch; either a full genotyping, or just a pileup depending on
    'args.pileup_only'; the results are written to 'filename'.
    """
    if args.pileup_only:
        return setup_mpileup_batch(args, regions, filename)
    return setup_genotyping_batch(args, regions, filename)


def poll_processes(procs, finalize=False):
    """Polls a set of processes, raising an exception if any have failed; if
    'finalize' is true, the function waits for each process to terminate.
    """
    for (key, proc) in procs.iteritems():
        returncode = (proc.wait() if finalize else proc.poll())
        if (finalize and returncode) or not (finalize or returncode is None):
            # Note that the proc.call property is set in function 'popen' above
            message = "Error with command %r, return-code = %s\n" \
                      "    Call = %r" % (key, returncode, proc.call)
            raise BatchError(message)


def run_batch((args, regions, filename)):
    setup = setup_batch(args, regions, filename)
    try:
        bam_handle_in = setup["handles"]["bam_in"]
        bam_handle_out = setup["handles"]["bam_out"]
        poll_processes(setup["procs"])

        nrecords = 0
        regions.reverse()
        while regions:
            region_aend = 0
            contig, start, end = regions[-1]
            for record in bam_handle_in.fetch(contig, start):
                current_aend = record.aend
                region_aend = max(region_aend, current_aend)
                if record.pos > end:
                    last_contig, _, _ = regions.pop()
                    if not regions:
                        break

                    contig, start, end = regions[-1]
                    if (region_aend + _BAM_BLOCK_SIZE < start) \
                            or (contig != last_contig):
                        break

                if current_aend >= start:
                    bam_handle_out.write(record)
                    if nrecords >= 100000:
                        poll_processes(setup["procs"])
                        nrecords = 0
                    nrecords += 1
            else:  # Reached the end of this contig
                while regions and (regions[-1][0] == contig):
                    regions.pop()

        setup["handles"]["bam_out"].close()

        poll_processes(setup["procs"], finalize=True)
        return filename
    except:
        traceback.print_exc()
        raise
    finally:
        cleanup_batch(setup)


###############################################################################
###############################################################################

def init_worker_thread():
    """Init function for subprocesses created by multiprocessing.Pool: Ensures
    that KeyboardInterrupts only occur in the main process, allowing us to do
    proper cleanup.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


###############################################################################
###############################################################################

def merge_bed_regions(regions):
    """Takes a sequence of bed regions [(contig, start, end), ...], which is
    assumed to be sorted by contig and coordiates, and returns a list in which
    overlapping records are merged into one larger region.
    """
    merged = []
    last_contig = last_start = last_end = None
    for record in regions:
        if (record.contig != last_contig) or (record.start > last_end):
            if last_contig is not None:
                merged.append((last_contig, last_start, last_end))
            last_contig = record.contig
            last_start = record.start
            last_end = record.end
        else:
            last_start = min(last_start or 0, record.start)
            last_end = max(last_end, record.end)

    if last_contig is not None:
        merged.append((last_contig, last_start, last_end))

    return merged


def create_batches(args, regions):
    """Yields a sequence of batches that may be passed to the 'run_batch'
    function; each batch consists of the 'args' object, a set of BED regions,
    and a destination filename. The set of BED regions is derived by splitting
    the total set of regions into args.nbatches portions.
    """
    tmpl = "{0}.batch_%03i".format(args.destination)
    def _get_batch_fname(count):
        """Returns a filename for batch number 'count'."""
        if count:
            return tmpl % (count,)
        return args.destination

    total_size = sum(end - start for (_, start, end) in regions)
    batch_size = total_size // args.nbatches + 5

    batch_count = 0
    current_batch = []
    current_total = 0
    for (contig, start, end) in regions:
        while (end - start) + current_total > batch_size:
            new_end = start + batch_size - current_total
            current_batch.append((contig, start, new_end))
            start = new_end
            yield args, current_batch, _get_batch_fname(batch_count)
            current_batch = []
            current_total = 0
            batch_count += 1
        current_batch.append((contig, start, end))
        current_total += end - start

    if current_batch:
        yield args, current_batch, _get_batch_fname(batch_count)


def merge_batch_results(filenames_iter):
    """Takes a multiprocessing.imap iterator yielding filenames of completed
    batches (gzipped vcf or mpileup files), and writes these into the
    file-handle out.
    """
    while True:
        try:
            # A timeout allows iteruption by the user, which is not the
            # case otherwise. The value is arbitrary.
            target_filename = filenames_iter.next(60)
            print "Merging into file: %r" % (target_filename,)
            break
        except multiprocessing.TimeoutError:
            pass
        except StopIteration:
            return

    with open(target_filename, "r+") as target_handle:
        while True:
            try:
                filename = filenames_iter.next(60)
                print "   - Processing batch: %r" % (filename,)

                # BGZip is terminated by 28b empty block (cf. ref)
                # While the standard implies that these should be ignored
                # if not actually at the end of the file, the tabix tool
                # stops processing at the first such block it encounters
                target_handle.seek(-28, 2)
                with open(filename) as input_handle:
                    shutil.copyfileobj(input_handle, target_handle)
                os.remove(filename)
            except multiprocessing.TimeoutError:
                pass
            except StopIteration:
                break


def collect_regions(args, bam_input_handle):
    """Returns the regions to be genotyped / pileup'd, as a list of bed-regions
    in the form (contig, start, end), where start is zero-based, and end is
    open based.
    """
    if args.bedfile is not None:
        regions = list(read_bed_file(args.bedfile))
        sort_bed_by_bamfile(bam_input_handle, regions)
        regions = merge_bed_regions(regions)
    else:
        regions = []
        for (name, length) in zip(bam_input_handle.references,
                                  bam_input_handle.lengths):
            regions.append((name, 0, length))
    return regions


def process_batches(args, batches):
    """Runs a set of batches, and merges the resulting output files if more
    than one batch is included.
    """
    nbatches = min(args.nbatches, len(batches))
    pool = multiprocessing.Pool(nbatches, init_worker_thread)

    try:
        batches = pool.imap(run_batch, batches, 1)
        merge_batch_results(batches)
    except:
        pool.terminate()
        pool.join()
        raise


def create_empty_bgz(destination):
    """Writes an empty BGZip file to the given destination; this file contains
    a single empty BGZip block (28b).
    """
    with open(destination, "w") as output:
        # Empty BGZip block
        output.write("\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42")
        output.write("\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00")
        output.write("\x00\x00")


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("bamfile", metavar='INPUT',
                        help="Sorted and indexed BAM file.")
    parser.add_argument("destination", metavar='OUTPUT',
                        help="BGZip compressed VCF or pileup. Also used as "
                             "prefix for temporary files.")
    parser.add_argument('--bedfile', default=None,
                        help="Optional bedfile, specifying regions to pileup "
                             "or genotype [Default: %(default)s].")
    parser.add_argument('--mpileup-argument', default=[], action="append",
                        help="Pass argument to 'samtools mpileup'; must be "
                             "used as follows: --mpileup-argument=-argument "
                             "for arguments without values, and "
                             "--mpileup-argument=-argument=value for "
                             "arguments with values.")
    parser.add_argument('--bcftools-argument', default=[], action="append",
                        help="Pass argument to 'bcftools view'; see the "
                             "--mpileup-argument command description.")
    parser.add_argument('--pileup-only', default=False, action="store_true",
                        help="Only run 'samtools mpileup', generating a text "
                             "pileup instead of a VCF file [Default: off].")
    parser.add_argument('--nbatches', metavar="N", default=1, type=int,
                        help="Split the BED into N number of batches, which "
                             "are run in parallel [Default: %(default)s].")
    parser.add_argument('--overwrite', default=False, action="store_true",
                        help="Overwrite output if it already exists "
                             "[Default: no].")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    if os.path.exists(args.destination) and not args.overwrite:
        sys.stderr.write("Output already exists; use --overwrite to allow "
                         "overwriting of this file.\n")
        return 1

    with pysam.Samfile(args.bamfile) as bam_input_handle:
        regions = collect_regions(args, bam_input_handle)
    batches = list(create_batches(args, regions))
    if not batches:
        create_empty_bgz(args.destination)
        return 0

    try:
        return process_batches(args, batches)
    except BatchError, error:
        sys.stderr.write("ERROR while processing BAM:\n")
        sys.stderr.write("    %s\n" % (error,))
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
