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
import pypeline
import pypeline.common.versions as versions

from pypeline.atomiccmd.builder import \
    AtomicCmdBuilder


def new(command, *args, **kwargs):
    return _COMMANDS[command](*args, **kwargs)


_version_eq = versions.EQ(*pypeline.__version_info__)
VERSION_PALEOMIX = versions.Requirement(call=["paleomix", "help"],
                                        search=r"v(\d+)\.(\d+)\.(\d+)",
                                        checks=_version_eq,
                                        priority=100)


def _build_cat_command():
    """Returns a AtomicCmdBuilder for the 'paleomix cat' command."""
    return AtomicCmdBuilder(["paleomix", "cat"],
                            EXEC_GZIP="gzip",
                            EXEC_BZIP="bzip2",
                            EXEC_CAT="cat",
                            CHECK_PALEOMIX=VERSION_PALEOMIX)


def _build_duphist_command():
    """Returns a AtomicCmdBuilder for the 'paleomix duphist' command."""
    return AtomicCmdBuilder(["paleomix", "duphist"],
                            CHECK_PALEOMIX=VERSION_PALEOMIX)


_COMMANDS = {
    "cat": _build_cat_command,
    "duphist": _build_duphist_command,
}
