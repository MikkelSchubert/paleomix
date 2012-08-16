#!/usr/bin/python
import os
import glob

from pypeline.common.fileutils import add_postfix

def prefix_to_filename(bwa_prefix):
    return os.path.splitext(os.path.basename(bwa_prefix))[0]


def sample_path(record, bwa_prefix = ""):
    """Returns """

    if bwa_prefix:
        bwa_prefix = "." + prefix_to_filename(bwa_prefix)
    return os.path.join("results", record["Name"] + bwa_prefix, record["Sample"])


def library_path(record, bwa_prefix = ""):
    """Returns """
    return os.path.join(sample_path(record, bwa_prefix), 
                        record["Library"])

def full_path(record, bwa_prefix = ""):
    """Returns """
    return os.path.join(library_path(record, bwa_prefix), 
                        record["Barcode"])



def target_path(record, bwa_prefix):
    filename = "%s.%s.bam" % (record["Sample"], prefix_to_filename(bwa_prefix))
    return os.path.join("results", filename)


def collect_files(record):
    """ """
    template = record["Path"]
    if template.format(Pair = 1) == template:
        return {"R1" : list(sorted(glob.glob(template))),
                "R2" : []}
    
    files = {}
    for (ii, name) in enumerate(("R1", "R2"), start = 1):
        files[name] = list(sorted(glob.glob(template.format(Pair = ii))))
    return files
