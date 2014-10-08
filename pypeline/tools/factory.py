#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 Mikkel Schubert <MSchubert@snm.ku.dk>
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
"""Factory for AtomicCmdBuilders for the various PALEOMIX commands.

Ensures
"""
import sys
import os.path
import pypeline
import pypeline.common.versions as versions

from pypeline.atomiccmd.builder import \
    AtomicCmdBuilder


# Get actual path of the 'paleomix' script used to invoke the current session;
# this is done to avoid conflict with local installs vs global installs.
_PALEOMIX_PATH = os.path.realpath(sys.argv[0])


def new(command, *args, **kwargs):
    """Returns AtomicCmdBuilder setup to call the tools accessible through the
    'paleomix' command-line tool. This builder adds executable / version checks
    for the specified command, but does not add any arguments. Thus, calling
    new with the argument "cat" produces the equivalent of ["paleomix", "cat"].
    """
    if command in _SPECIAL_COMMANDS:
        return _SPECIAL_COMMANDS[command](*args, **kwargs)
    return _build_generic_command(command)


_VERSION_EQ = versions.EQ(*pypeline.__version_info__)
VERSION_PALEOMIX = versions.Requirement(call=[_PALEOMIX_PATH, "help"],
                                        search=r"v(\d+)\.(\d+)\.(\d+)",
                                        checks=_VERSION_EQ,
                                        priority=100)


def _build_generic_command(argument):
    """Returns a AtomicCmdBuilder for a regular 'paleomix ...' command."""
    return AtomicCmdBuilder([_PALEOMIX_PATH, argument],
                            CHECK_PALEOMIX=VERSION_PALEOMIX)


def _build_cat_command():
    """Returns a AtomicCmdBuilder for the 'paleomix cat' command."""
    return AtomicCmdBuilder([_PALEOMIX_PATH, "cat"],
                            EXEC_GZIP="gzip",
                            EXEC_BZIP="bzip2",
                            EXEC_CAT="cat",
                            CHECK_PALEOMIX=VERSION_PALEOMIX)


def _build_create_pileup_command(outfile):
    """Returns a AtomicCmdBuilder for a regular 'paleomix ...' command."""
    return AtomicCmdBuilder([_PALEOMIX_PATH, "create_pileup", outfile],
                            CHECK_PALEOMIX=VERSION_PALEOMIX)


_SPECIAL_COMMANDS = {
    "cat": _build_cat_command,
    "create_pileup": _build_create_pileup_command,
}
