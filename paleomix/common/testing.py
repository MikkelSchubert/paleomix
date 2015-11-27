#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import os
import sys
import shutil
import tempfile

import nose
from nose.tools import assert_equal

from paleomix.common.fileutils import \
     make_dirs


def assert_list_equal(iter_a, iter_b):
    """Compare two values, after first converting them to lists.
    This ensures that lazily generated results can be compared."""
    list_a = list(iter_a)
    list_b = list(iter_b)

    assert_equal(list_a, list_b)


def with_temp_folder(func):
    """Decorator for unit-tests:
    Creates a unique temporary folder before running 'func'. The
    function is is assumed to take at least one parameter, the first
    of which is assumed to represent the temporary folder."""
    temp_root = os.path.join(tempfile.gettempdir(), os.getlogin())
    make_dirs(temp_root) # Ensure that this subdirectory exists

    @nose.tools.istest
    def _wrapper(*args, **kwargs):
        try:
            temp_folder = None
            temp_folder = tempfile.mkdtemp(dir    = temp_root,
                                           prefix = "paleomix_unit")
            func(temp_folder, *args, **kwargs)
        finally:
            if temp_folder:
                shutil.rmtree(temp_folder)
    _wrapper.__name__ = func.__name__ + "__wrapped_by_with_temp_folder"
    return _wrapper


class Monkeypatch:
    """Replaces a function/object in a module with the specified wrapper
     upon entry, reverting the change upon exit from the with statement.
     A full path to the given function is required, for example
       'os.path.join'."""
    def __init__(self, path, wrapper):
        self.wrapper = wrapper

        parts = path.split(".")
        assert len(parts) > 1
        self.module, self.object  = None, sys.modules[parts[0]]
        for path_cmp in parts[1:]:
            self.module, self.object = self.object, getattr(self.object, path_cmp)
        self.name = parts[-1]

    def __enter__(self):
        setattr(self.module, self.name, self.wrapper)
        return self

    def __exit__(self, _type, _value, _traceback):
        setattr(self.module, self.name, self.object)


class SetWorkingDirectory:
    """Sets the current working directory upon entry to that specified,
    in the constructor upon entry, and reverts to the previously used
    directory upon exiting a with statement."""
    def __init__(self, path):
        self._old_cwd = None
        self._new_cwd = path

    def __enter__(self):
        self._old_cwd = os.getcwd()
        os.chdir(self._new_cwd)

    def __exit__(self, _type, _value, _traceback):
        os.chdir(self._old_cwd)


def set_file_contents(fname, contents):
    with open(fname, "w") as handle:
        handle.write(contents)


def get_file_contents(fname):
    with open(fname) as handle:
        return handle.read()
