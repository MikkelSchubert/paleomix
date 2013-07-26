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
import time
import shutil
import atexit

import nose
from nose.tools import with_setup, assert_equal


# Simple support for assert_in in Python v2.6
try:
    from nose.tools import assert_in # pylint: disable=W0611
except ImportError:
    def assert_in(item, lst):
        assert item in lst

def assert_list_equal(iter_a, iter_b):
    """Compare two values, after first converting them to lists.
    This ensures that lazily generated results can be compared."""
    list_a = list(iter_a)
    list_b = list(iter_b)

    assert_equal(list_a, list_b)




EPOCH = time.time()
TEST_ROOT = "/tmp/%s/pypeline/tests_%i" % (os.getlogin(), EPOCH,)
UNIT_TEST_ROOT = os.path.join(TEST_ROOT, "unit")
FUNC_TEST_ROOT = os.path.join(TEST_ROOT, "func")

_COUNTER = 0
def with_temp_folder(func):
    """Creates a unique temporary folder before running 'func'. The
    function is is assumed to take at least one parameter, the first
    of which is assumed to represent the temporary folder."""
    global _COUNTER # pylint: disable=W0603
    tmp_root = os.path.join(UNIT_TEST_ROOT, "%03i_%s" % (_COUNTER, func.__name__,))
    _COUNTER += 1

    def _setup():
        os.makedirs(tmp_root)

    def _teardown():
        shutil.rmtree(tmp_root)

    @nose.tools.istest
    def _wrapper(*args, **kwargs):
        return func(tmp_root, *args, **kwargs)
    _wrapper.__name__ = func.__name__ + "__wrapped_by_with_temp_folder"

    return with_setup(_setup, _teardown)(_wrapper)



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


class RequiredCall(Monkeypatch):
    """Monkey-patches the specified function, and checks that this function has
    been called upon exiting from a with statement. If 'args' or 'kwargs' are
    specified, the caller must match these in order to be considered a valid
    call to the original function. A full path to the given function is
    required, for example: 'os.path.join'."""
    def __init__(self, path, args = None, kwargs = None):
        self._path = path
        self._args = tuple(args) if args else args
        self._kwargs = dict(kwargs) if kwargs else kwargs
        self._was_called = False
        Monkeypatch.__init__(self, path, self._called)

    def _called(self, *args, **kwargs):
        if (self._args is None) or (self._args == args):
            if (self._kwargs is None) or (self._kwargs == kwargs):
                self._was_called = True
        return self.object(*args, **kwargs)

    def __exit__(self, type_, value, traceback):
        Monkeypatch.__exit__(self, type_, value, traceback)
        if not value:
            assert self._was_called, "%s was not called with required arguments!" % self._path


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



@atexit.register
def _cleanup_temp_folders():
    # Avoid (when possible) leaving empty folders behind after each run
    for root in (FUNC_TEST_ROOT, UNIT_TEST_ROOT, TEST_ROOT):
        try:
            os.rmdir(root)
        except OSError:
            pass
