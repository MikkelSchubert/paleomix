#!/usr/bin/python
import os
import sys

from pypeline.common.text import padded_table
                     
def read_alignment_records(root):
    with open(os.path.join(root, "SampleSheet.csv")) as records:
        header = records.readline().strip().split(",")
        for line in records:
            yield dict(zip(header, line.strip().split(",")))


def main(argv):
    lines = []
    for root in argv:
        for record in read_alignment_records(root):
            record["Lane"] = int(record["Lane"])
            line = [record["SampleID"],
                    record["SampleID"],
                    record["Index"],
                    record["FCID"],
                    "Illumina",
                    os.path.join(root, "%(SampleID)s_%(Index)s_L%(Lane)03i_R{Pair}_*.fastq.gz" % record)]
            lines.append(line)
    lines.sort(key = lambda v: map(str.lower, v))
    lines.insert(0, ["Name", "Sample", "Library", "Barcode", "Platform", "Path"])

    for line in padded_table(lines):
        print line


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
