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
"""


"""
import os
import sys
import gzip
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


class BatchError(RuntimeError):
    pass


# Size of smallest block in (linear) BAM index (= 2 << 14)
_BAM_BLOCK_SIZE = 16384


###############################################################################
###############################################################################

def popen(call, *args, **kwargs):
    """Equivalent to """
    proc = subprocess.Popen(call, *args, **kwargs)
    proc.call = tuple(call)
    return proc


###############################################################################
###############################################################################
## CLI functions

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
## Batch functions
##

def can_reduce_batch(args, regions):
    first = regions[0]
    last = regions[-1]

    if first[0] != last[0]:  # Same contig?
        return False

    span = last[-1] - first[1]
    nbases = sum((end - start) for (_, start, end) in regions)

    return 1.0 - nbases / float(span) < args.max_gappiness


def reduce_batch(args, regions):
    if can_reduce_batch(args, regions):
        first = regions[0]
        last = regions[-1]

        return "%s:%i-%i" % (first[0], first[1], last[-1])


###############################################################################
###############################################################################
## Common functions

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


def write_sample_file(prefix, samples):
    fpath = prefix + ".samples"
    with open(fpath, "w") as sample_handle:
        for sample in sorted(samples):
            sample_handle.write("%s     1\n" % (sample,))
    return fpath


def make_bam_pipe(prefix):
    fpath_pipe = prefix + ".pipe"
    if os.path.exists(fpath_pipe):
        os.remove(fpath_pipe)
    os.mkfifo(fpath_pipe)
    return fpath_pipe


def setup_basic_batch(args, regions, prefix, bgzip, func):
    setup = {"files": {},
             "temp_files": {},
             "procs": {},
             "handles": {}}

    try:
        reduced_region = reduce_batch(args, regions)
        setup["files"]["bed"] = write_bed_file(prefix, regions)
        setup["temp_files"]["bed"] = setup["files"]["bed"]

        if reduced_region:
            print "Processing batch in region %r" % (reduced_region,)
            setup["files"]["pipe"] = args.bamfile
        else:
            setup["files"]["pipe"] = make_bam_pipe(prefix)
            setup["temp_files"]["pipe"] = setup["files"]["pipe"]

        stdout = func(setup, reduced_region)

        # Quick compression if intermediate result, bgzip otherwise
        gzip_call = ["bgzip"] if bgzip else ["gzip", "-1"]
        setup["handles"]["outfile"] = open(prefix, "w")
        zip_proc = popen(gzip_call,
                         stdin=stdout,
                         stdout=setup["handles"]["outfile"],
                         close_fds=True)

        setup["procs"]["gzip"] = zip_proc

        if not reduced_region:
            setup["handles"]["bam_in"] = pysam.Samfile(args.bamfile)
            setup["handles"]["bam_out"] = \
                pysam.Samfile(setup["files"]["pipe"], "wbu",
                              template=setup["handles"]["bam_in"])

        if reduced_region:
            setup["files"].pop("pipe")

        return setup
    except:
        traceback.print_exc()
        cleanup_batch(setup)
        raise


###############################################################################
###############################################################################
## Pileup batch generation

def setup_mpileup_batch(args, regions, prefix, bgzip):
    def _create_mpileup_proc(setup, region):
        mpileup_args = {"-l": setup["files"]["bed"]}
        if region is not None:
            mpileup_args["-r"] = region

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

    return setup_basic_batch(args, regions, prefix, bgzip,
                             _create_mpileup_proc)


###############################################################################
###############################################################################
## Genotyping batch generation

def setup_genotyping_batch(args, regions, prefix, bgzip):
    def _create_genotyping_proc(setup, region):
        mpileup_args = {"-u": None,
                        "-l": setup["files"]["bed"]}
        if region is not None:
            mpileup_args["-r"] = region

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

    return setup_basic_batch(args, regions, prefix, bgzip,
                             _create_genotyping_proc)


###############################################################################
###############################################################################

def setup_batch(args, regions, filename, bgzip):
    if args.pileup_only:
        return setup_mpileup_batch(args, regions, filename, bgzip)
    return setup_genotyping_batch(args, regions, filename, bgzip)


def poll_processes(procs, finalize=False):
    for (key, proc) in procs.iteritems():
        returncode = (proc.wait() if finalize else proc.poll())
        if (finalize and returncode) or not (finalize or returncode is None):
            message = "Error with command %r, return-code = %s\n" \
                      "    Call = %r" % (key, returncode, proc.call)
            raise BatchError(message)


def run_batch((args, regions, filename), bgzip=False):
    setup = setup_batch(args, regions, filename, bgzip)
    try:
        if "bam_out" in setup["handles"]:
            bam_handle_in = setup["handles"]["bam_in"]
            bam_handle_out = setup["handles"]["bam_out"]

            nrecords = 0
            regions.reverse()
            while regions:
                #aend = 0
                contig, start, end = regions[-1]
                for record in bam_handle_in.fetch(contig, start):
                    bam_handle_out.write(record)
                    if nrecords >= 100000:
                        poll_processes(setup["procs"])
                        nrecords = 0
                    nrecords += 1

                    #aend = max(aend, record.aend)
                    if record.pos > end:
                        last_contig, _, _ = regions.pop()
                        if not regions:
                            break

                        contig, start, end = regions[-1]
                        # Only do a new fetch if the next region is
                        # likely to be located in a different block
                        if (record.pos + _BAM_BLOCK_SIZE < start) \
                                or (contig != last_contig):
                            break
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
    proper cleanup."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


###############################################################################
###############################################################################

def batch_size(regions):
    return sum((end - start) for (_, start, end) in regions)


def merge_bed_regions(regions):
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
    merged.append((last_contig, last_start, last_end))
    return merged


def create_batches(args, regions):
    tmpl = "{0}.batch_%03i".format(args.destination)

    total_size = sum(end - start for (_, start, end) in regions)
    batch_size = total_size // args.nbatches + 5

    batch_count = 0
    current_batch = []
    current_total = 0
    for (contig, start, end) in regions:
        while ((end - start) + current_total > batch_size):
            new_end = start + batch_size - current_total
            current_batch.append((contig, start, new_end))
            start = new_end
            yield args, current_batch, tmpl % batch_count
            current_batch = []
            current_total = 0
            batch_count += 1
        current_batch.append((contig, start, end))
        current_total += end - start

    if current_batch:
        yield args, current_batch, tmpl % batch_count


def merge_batch_results(filenames_iter, proc):
    has_header = False
    while True:
        try:
            # A timeout allows iteruption by the user, which is not the
            # case otherwise. The value is arbitrary.
            filename = filenames_iter.next(60)
            print "Processing batch %r" % (filename,)
            with gzip.open(filename) as input_handle:
                if has_header:  # Skip header in all but the first file
                    line = input_handle.readline()
                    while line and line.startswith("#"):
                        line = input_handle.readline()
                    proc.stdin.write(line)

                shutil.copyfileobj(input_handle, proc.stdin)
                has_header = True
            os.remove(filename)
        except multiprocessing.TimeoutError:
            pass
        except StopIteration:
            break


def collect_regions(args, bam_input_handle):
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
    nbatches = min(args.nbatches, len(batches))
    pool = multiprocessing.Pool(nbatches, init_worker_thread)

    proc = None
    with open(args.destination, "w") as output_handle:
        try:
            proc = popen(["bgzip"],
                         stdin=subprocess.PIPE,
                         stdout=output_handle,
                         close_fds=True)

            batches = pool.imap(run_batch, batches, 1)
            merge_batch_results(batches, proc)

            proc.stdin.close()
            return proc.wait()
        except:
            pool.terminate()
            pool.join()
            proc.terminate()
            proc.wait()
            raise


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("bamfile", metavar='INPUT')
    parser.add_argument("destination", metavar='OUTPUT')
    parser.add_argument('--mpileup-argument', default=[], action="append")
    parser.add_argument('--bcftools-argument', default=[], action="append")
    parser.add_argument('--bedfile', default=None)
    parser.add_argument('--max-gappiness', default=0.2, type=float)
    parser.add_argument('--pileup-only', default=False, action="store_true")
    parser.add_argument('--nbatches', default=1, type=int)
    # Max gappiness
    # Per contig

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    with pysam.Samfile(args.bamfile) as bam_input_handle:
        regions = collect_regions(args, bam_input_handle)
    batches = list(create_batches(args, regions))

    try:
        if len(batches) > 1:
            return process_batches(args, batches)
        else:
            batch = list(batches[0])
            batch[-1] = args.destination
            if not run_batch(batch, bgzip=True):
                return 1
    except BatchError, error:
        sys.stderr.write("ERROR while processing BAM:\n")
        sys.stderr.write("    %s\n" % (error,))
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
