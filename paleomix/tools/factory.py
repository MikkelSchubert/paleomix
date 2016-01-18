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

Ensures that the version called corresponds to the running version, in case
multiple versions are present in the users' PATH, or that the current version
is not available from the users' PATH.
"""
import sys

import paleomix.main

from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder


def new(command, *args, **kwargs):
    """Returns AtomicCmdBuilder setup to call the tools accessible through the
    'paleomix' command-line tool. This builder adds executable / version checks
    for the specified command, but does not add any arguments. Thus, calling
    new with the argument "cat" produces the equivalent of ["paleomix", "cat"].
    """
    if command in _SPECIAL_COMMANDS:
        return _SPECIAL_COMMANDS[command](*args, **kwargs)
    return _build_paleomix_command(command, *args, **kwargs)


def _build_cat_command():
    """Returns an AtomicCmdBuilder for the 'paleomix cat' command."""
    return _build_paleomix_command("cat",
                                   EXEC_GZIP="gzip",
                                   EXEC_BZIP="bzip2",
                                   EXEC_CAT="cat")


def _build_paleomix_command(*args, **kwargs):
    """Returns an AtomicCmdBuilder for a regular 'paleomix ...' command."""
    interpreter = sys.executable
    script = paleomix.main.__file__

    return AtomicCmdBuilder((interpreter, script) + args,
                            AUX_PALEOMIX=script,
                            **kwargs)


_SPECIAL_COMMANDS = {
    "cat": _build_cat_command,
}
