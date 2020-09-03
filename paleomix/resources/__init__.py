#!/usr/bin/python3
#
# Copyright (c) 2016 Mikkel Schubert <MikkelSch@gmail.com>
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
import logging
import os
import shutil

from pkg_resources import Requirement, cleanup_resources, resource_filename


_REQUIREMENT = Requirement.parse("PALEOMIX")


def rscript(tool, script):
    """Returns the path to an Rscript for a given tool."""

    path = os.path.join("paleomix", "resources", "rscripts", tool, script)

    return resource_filename(_REQUIREMENT, path)


def report(tool, filename):
    """Returns the path to a report-file for a given tool."""

    path = os.path.join("paleomix", "resources", "reports", tool, filename)

    return resource_filename(_REQUIREMENT, path)


def add_copy_example_command(subparsers):
    parser = subparsers.add_parser("example", help="Create example project")

    parser.add_argument(
        "destination",
        default=".",
        nargs="?",
        help="Destination folder for example data.",
    )


def copy_example(tool, args):
    """Command-line interface to copy a folder containing example data to a
    folder specified by the user. Arguments are a tool name (e.g.
    'bam_pipeline'), and any command-line options specified by the user;
    returns 0 on success, or 1 on errors.
    """
    log = logging.getLogger(__name__)
    destination = os.path.join(args.destination, tool)
    log.info("Copying example project to %r", args.destination)

    if os.path.exists(destination):
        log.error("Example folder already exists at %r", destination)
        return 1

    try:
        source = resource_filename("paleomix.resources", os.path.join("examples", tool))

        shutil.copytree(source, destination)
    finally:
        cleanup_resources()

    log.info("Sucessfully saved example in %r", destination)

    return 0
