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
import subprocess

import pypeline.ui
import pypeline.atomiccmd

from pypeline.common.utilities import safe_coerce_to_tuple


def select_cat(filename):
     with open(filename) as source:
          header = source.read(2)
          if header == "\x1f\x8b":
               return "zcat"
          elif header == "BZ":
               return "bzcat"

          return "cat"




def main(input_files, output_file):
     executable = select_cat(input_files[0])
     command = [executable] + input_files
     
     with open(output_file, "w") as output:
          proc = subprocess.Popen(command, 
                                  stdout    = output, 
                                  close_fds = True)
          return proc.wait()


def build_atomiccmd(cls, input_files, output_file):
     filename = os.path.realpath(__file__)
     if filename.endswith(".pyc"):
          filename = filename[:-1]

     command  = [filename, "--output", "%(TEMP_OUT_FILE)s"]
     files = {"TEMP_OUT_FILE" : output_file}
     for (ii, filename) in enumerate(safe_coerce_to_tuple(input_files), start = 1):
          files["IN_FILE_%02i" % ii] = filename
          command.append("%%(IN_FILE_%02i)s" % ii)
          
     return cls(command, **files)
     

if __name__ == '__main__':
     parser = optparse.OptionParser()
     parser.add_option("--output")

     config, args = parser.parse_args()
     if not args:
          pypeline.ui.print_err("%s takes at least one input file" % sys.argv[0])
          sys.exit(1)
     elif not config.output:
          pypeline.ui.print_err("%s required that --output is set" % sys.argv[0])
          sys.exit(1)
                               
     main(input_files   = args,
          output_file   = config.output)
