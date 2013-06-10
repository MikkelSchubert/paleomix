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

from nose.tools import with_setup


EPOCH = time.time()
TEST_ROOT = "tests/runs/unit_%i/temp/" % (EPOCH,)


def with_temp_folder(func):
    """Creates a unique temporary folder before running 'func'. The
    function is is assumed to take a single parameter, representing
    the path to the newly created folder"""
    tmp_root = os.path.join(TEST_ROOT, func.__name__)
    def _setup():
        os.makedirs(tmp_root)

    def _teardown():
        shutil.rmtree(tmp_root)

    def _wrapper():
        return func(tmp_root)
    _wrapper.__name__ = func.__name__ + "__wrapped_by_with_temp_folder"
    return with_setup(_setup, _teardown)(_wrapper)



class monkeypatch:
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


class require_call(monkeypatch):
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
        monkeypatch.__init__(self, path, self._called)

    def _called(self, *args, **kwargs):
        if (self._args is None) or (self._args == args):
            if (self._kwargs is None) or (self._kwargs == kwargs):
                self._was_called = True
        return self.object(*args, **kwargs)

    def __exit__(self, type, value, traceback):
        monkeypatch.__exit__(self, type, value, traceback)
        if not value:
            assert self._was_called, "%s was not called with required arguments!" % self._path


class set_cwd:
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
