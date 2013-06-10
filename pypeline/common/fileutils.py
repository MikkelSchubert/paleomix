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
import uuid
import shutil
import errno
import bz2
import gzip

from pypeline.common.utilities import safe_coerce_to_tuple


def add_postfix(filename, postfix):
    """Ads a postfix to a filename (before any extensions that filename may have)."""
    filename, ext = os.path.splitext(filename)
    return filename + postfix + ext


def swap_ext(filename, ext):
    """Replaces the existing extension of a filename with the specified extension,
    ensuring that the extension is prefixed by a '.'. If no extension is specified,
    other than potentially a dot, the existing extension is stripped."""
    filename, _ = os.path.splitext(filename)
    if ext in ("", "."):
        return filename

    if not ext.startswith("."):
        ext = "." + ext

    return filename + ext


def reroot_path(root, filename):
    """Returns the basename of filename, joined to root."""
    return os.path.join(root, os.path.basename(filename))


def create_temp_dir(root):
    """Creates a temporary directory, accessible only by the owner,
    at the specified location. The folder name is randomly generated,
    and only the current user has access"""
    def _generate_path():
        return os.path.join(root, str(uuid.uuid4()))

    path = _generate_path()
    while not make_dirs(path, mode = 0700):
        path = _generate_path()
    return path


def missing_files(filenames):
    """Given a list of filenames, returns a list of those that
    does not exist. Note that this function does not differentiate
    between files and folders."""
    result = []
    for filename in safe_coerce_to_tuple(filenames):
        if not os.path.exists(filename):
            result.append(filename)

    return result


def modified_after(younger, older):
    """Returns true any of the files expected to be 'younger' have
    been modified after any of the files expected to be 'older'."""
    def get_mtimes(filenames):
        for filename in filenames:
            yield os.path.getmtime(os.path.realpath(filename))

    younger_time = max(get_mtimes(safe_coerce_to_tuple(younger)))
    older_time   = min(get_mtimes(safe_coerce_to_tuple(older)))

    return younger_time > older_time


def is_executable(filename):
    """Returns true if the specified file is an executable file."""
    return os.path.isfile(filename) and os.access(filename, os.X_OK)


def executable_exists(filename):
    """Returns true if the filename refers to an executable file,
    either by relative or full path, or if the executable is found
    on the current PATH."""
    if os.path.dirname(filename):
        return is_executable(filename)

    for path in os.environ["PATH"].split(os.pathsep):
        if is_executable(os.path.join(path, filename)):
            return True

    return False


def missing_executables(filenames):
    result = []
    for filename in safe_coerce_to_tuple(filenames):
        if not executable_exists(filename):
            result.append(filename)
    return result


def make_dirs(directory, mode = 0777):
    """Wrapper around os.makedirs to make it suitable for using
    in a multithreaded/multiprocessing enviroment: Unlike the
    regular function, this wrapper does not throw an exception if
    the directory already exists, which may happen if another
    thread/process created the directory during the function call.

    Returns true if a new directory was created, false if it
    already existed. Other errors result in exceptions."""
    if not directory:
        raise ValueError("Empty directory passed to make_dirs()")

    try:
        os.makedirs(directory, mode = mode)
        return True
    except OSError, error:
        # make_dirs be called by multiple subprocesses at the same time,
        # so only raise if the actual creation of the folder failed
        if error.errno != errno.EEXIST:
            raise
        return False


def move_file(source, destination):
    """Wrapper around shutils which ensures that the
    destination directory exists before moving the file."""
    dirname = os.path.dirname(destination)
    if dirname:
        make_dirs(dirname)
    shutil.move(source, destination)


def copy_file(source, destination):
    """Wrapper around shutils which ensures that the
    destination directory exists before copying the file."""
    dirname = os.path.dirname(destination)
    if dirname:
        make_dirs(dirname)
    shutil.copy(source, destination)



def open_ro(filename):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(filename)
    try:
        header = handle.read(2)
        handle.seek(0)

        if header == "\x1f\x8b":
            handle.close()
            # TODO: Re-use handle (fileobj)
            handle = gzip.open(filename)
        elif header == "BZ":
            handle.close()
            handle = bz2.BZ2File(filename)

        return handle
    except:
        handle.close()
        raise
