#!/usr/bin/python
"""
An alternative to 'samtools view -Sb -', needed because samtools does not
properly handle SAM-files with no records. For example, the following
command will not work:
$ samtools view -H INPUT | samtools view -Sb -
"""

import os
import sys
import optparse

import pysam

import pypeline.ui
import pypeline.atomiccmd


def main(min_quality, exclude_flags):
     input_file = pysam.Samfile("-", "r")
     output_file = pysam.Samfile("-", "wb", template = input_file)
     for read in input_file:
          if read.flag & exclude_flags:
               continue
          elif read.mapq < min_quality:
               continue
               
          output_file.write(read)
     output_file.close()
     input_file.close()


def build_atomiccmd(cls, min_quality, exclude_flags, stdin, output_file):
     filename = os.path.realpath(__file__)
     if filename.endswith(".pyc"):
          filename = filename[:-1]

     return cls([filename, 
                 "-q", min_quality,
                 "-F", exclude_flags],
                IN_STDIN   = stdin,
                OUT_STDOUT = output_file)
     

if __name__ == '__main__':
     parser = optparse.OptionParser()
     parser.add_option("-q", "--min-quality", type = int, default = 0)
     parser.add_option("-F", "--exclude-flags", type = int, default = 0)

     config, args = parser.parse_args()
     if args:
          pypeline.ui.print_err("%s does not take any arguments." % sys.argv[0])
          sys.exit(1)
                               
     main(min_quality   = config.min_quality,
          exclude_flags = config.exclude_flags)
