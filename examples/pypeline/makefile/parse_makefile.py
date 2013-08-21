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
"""Example usage of makefile support in the Pypeline libraries.

See documentation in pypeline.common.makefile for more information."""
import sys
import pprint

from pypeline.common.makefile import \
     MakefileError, \
     read_makefile, \
     WithoutDefaults, \
     Or, \
     IsInt, \
     IsFloat, \
     IsStr, \
     StringStartsWith


_SPECIFICATION_OF_OPTIONS = {
    StringStartsWith("--") : Or(IsInt, IsFloat),
    "--min-depth"          : IsInt(default = 8),
    "--max-depth"          : IsInt(default = 100),
}


_MAKEFILE_SPECIFICATION = {
    "Defaults"  : _SPECIFICATION_OF_OPTIONS,

    "VCF_Files" : {
        IsStr : {
            "Output_File" : IsStr,
            "Options" : WithoutDefaults(_SPECIFICATION_OF_OPTIONS),
        }
    }
}



def main(argv):
    for filename in argv:
        try:
            print "Reading makefile at %r" % (filename,)
            pprint.pprint(read_makefile(filename, _MAKEFILE_SPECIFICATION))
        except MakefileError, error:
            print "Error parsing makefile:\n%s" % (error,)
        print

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
