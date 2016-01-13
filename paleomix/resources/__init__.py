#!/usr/bin/python
#
# Copyright (c) 2016 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import argparse
import os
import shutil
import sys

from pkg_resources import Requirement, resource_filename


_REQUIREMENT = Requirement.parse("PALEOMIX")


def rscript(tool, script):
    """Returns the path to an Rscript for a given tool."""

    path = os.path.join("paleomix", "resources", "rscripts", tool, script)

    return resource_filename(_REQUIREMENT, path)


def report(tool, filename):
    """Returns the path to a report-file for a given tool."""

    path = os.path.join("paleomix", "resources", "reports", tool, filename)

    return resource_filename(_REQUIREMENT, path)


def copy_example(tool, argv):
    """Command-line interface to copy a folder containing example data to a
    folder specified by the user. Arguments are a tool name (e.g.
    'bam_pipeline'), and any command-line options specified by the user;
    returns 0 on success, or 1 on errors.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('root', help="Destination folder for example data.")

    args = parser.parse_args(argv)

    destination = os.path.join(args.root, tool)
    if os.path.exists(destination):
        sys.stderr.write("Example folder already exists at destination, "
                         "cannot proceed:\n")
        sys.stderr.write("  - %r\n" % (destination,))
        return 1

    path = os.path.join("paleomix", "resources", "examples", tool)
    source = resource_filename(_REQUIREMENT, path)

    shutil.copytree(source, destination)

    sys.stderr.write("Sucessfully saved example in %r\n" % (destination,))

    return 0
