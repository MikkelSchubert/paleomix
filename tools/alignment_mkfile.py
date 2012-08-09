#!/usr/bin/python
import os
import sys


_FIELD_SIZE = 7
def calc_field_size(field):
    return (len(field) + _FIELD_SIZE - 1) / _FIELD_SIZE


def calc_ntabs(lines):
    """Calculates the size in tabs for each column."""
    tabs = [0] * len(lines[0])
    for line in lines:
        for (ii, field) in enumerate(line):
            tabs[ii] = max(tabs[ii], calc_field_size(field))
    tabs[-1] = 0
    return tabs


def print_padded_line(tabs, line):
    """Pads a line, Using a list of sizes calculated by 'calc_ntabs'."""
    padded = []
    for (ii, field) in enumerate(line):
        padded.append(field + "\t" * (tabs[ii] - calc_field_size(field)))

    print "\t".join(padded)
                       


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
                    os.path.join(root, "%(SampleID)s_%(Index)s_L%(Lane)03i_%%(Pair)s_*.fastq.gz" % record)]
            lines.append(line)
    lines.sort(key = lambda v: map(str.lower, v))
    lines.insert(0, ["Name", "Sample", "Library", "Barcode", "Platform", "Path"])

    tabs = calc_ntabs(lines)
    for line in lines:
        print_padded_line(tabs, line)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
