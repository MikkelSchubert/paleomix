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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
#
import os
import glob

ROOT = "results"


def prefix_to_filename(bwa_prefix):
    return os.path.splitext(os.path.basename(bwa_prefix))[0]


def prefix_path(record, bwa_prefix = ""):
    bwa_prefix = prefix_to_filename(bwa_prefix)
    return os.path.join(ROOT, record["Name"], bwa_prefix)


def sample_path(record, bwa_prefix = ""):
    """Returns """

    return os.path.join(prefix_path(record, bwa_prefix), record["Sample"])


def library_path(record, bwa_prefix = ""):
    """Returns """
    return os.path.join(sample_path(record, bwa_prefix), 
                        record["Library"])

def full_path(record, bwa_prefix = ""):
    """Returns """
    return os.path.join(library_path(record, bwa_prefix), 
                        record["Barcode"])



def target_path(record, bwa_prefix):
    filename = "%s.%s.bam" % (record["Name"], prefix_to_filename(bwa_prefix))
    return os.path.join(ROOT, filename)

def mapdamage_path(name, library, bwa_prefix):
    dirname = "%s.%s.mapDamage" % (name, prefix_to_filename(bwa_prefix))
    return os.path.join(ROOT, dirname, library)



def is_paired_end(record):
    template = record["Path"]

    return (template.format(Pair = 1) != template)


def collect_files(record):
    """ """
    template = record["Path"]
    if not is_paired_end(record):
        return {"SE" : list(sorted(glob.glob(template)))}
    
    files = {}
    for (ii, name) in enumerate(("PE_1", "PE_2"), start = 1):
        files[name] = list(sorted(glob.glob(template.format(Pair = ii))))
    return files


def reference_sequence(bwa_prefix):
    if os.path.exists(bwa_prefix + ".fasta"):
        return bwa_prefix + ".fasta"
    elif os.path.exists(bwa_prefix + ".fa"):
        return bwa_prefix + ".fa"

    wo_ext = os.path.splitext(bwa_prefix)[0]
    if bwa_prefix != wo_ext:
        return reference_sequence(wo_ext)
