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
"""Common version requirements for pypeline/pipelines."""
import sys

import pypeline.common.versions as versions
from pypeline.common.console import print_err


def check_pypeline_requirements():
    for (module, url) in sorted(_REQUIRED_MODULES.items()):
        try:
            __import__(module)
        except ImportError, error:
            print_err("ERROR: Could not import required module '%s':\n\t- %s\n" % (module, error), file = sys.stderr)
            print_err("       If module is not installed, please download from '%s'.\n" % (url,), file = sys.stderr)
            print_err("       A local install may be performed using the following command:\n", file = sys.stderr)
            print_err("       $ python setup.py install --user\n\n", file = sys.stderr)
            return False

    for version in (_PYSAM_VERSION,):
        try:
            version()
        except versions.VersionRequirementError, error:
            print_err("Version requirements check failed for %s:\n   %s\n" \
                      % (version.name, error), file = sys.stderr)
            return False
    return True


# List of modules required to run the pipelines
_REQUIRED_MODULES = {"pysam" : "http://code.google.com/p/pysam/"}


# Cannot be a lambda due to need to be able to pickle function
def _get_pysam_version():
    return __import__("pysam").__version__
_PYSAM_VERSION = versions.Requirement(name   = "module 'pysam'",
                                     call   = _get_pysam_version,
                                     search = r"(\d+)\.(\d+)\.(\d+)",
                                     checks = versions.GE(0, 7, 5))
