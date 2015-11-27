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
import bz2
import gzip
import uuid
import errno
import types
import shutil

from paleomix.common.utilities import safe_coerce_to_tuple, \
     safe_coerce_to_frozenset


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
    while not make_dirs(path, mode = 0750):
        path = _generate_path()
    return path


def missing_files(filenames):
    """Given a list of filenames, returns a list of those that
    does not exist. Note that this function does not differentiate
    between files and folders."""
    result = []
    for filename in safe_coerce_to_frozenset(filenames):
        if not os.path.exists(filename):
            result.append(filename)

    return result


def modified_after(younger, older):
    """Returns true any of the files expected to be 'younger' have
    been modified after any of the files expected to be 'older'."""
    def get_mtimes(filenames):
        for filename in filenames:
            yield os.path.getmtime(filename)

    younger_time = max(get_mtimes(safe_coerce_to_frozenset(younger)))
    older_time   = min(get_mtimes(safe_coerce_to_frozenset(older)))

    return younger_time > older_time


def is_executable(filename):
    """Returns true if the specified file is an executable file."""
    return os.path.isfile(filename) and os.access(filename, os.X_OK)


def which_executable(filename):
    """Returns the path of the first executable in the PATH which
    matches the filename, or None if no match was found. If the
    filename contains a directory component, only that path is
    tested, and None is returned if that file is not an executable."""
    if os.path.dirname(filename):
        if is_executable(filename):
            return filename
        return None

    path_variable = os.environ.get("PATH")
    if not path_variable:
        return None

    for path in path_variable.split(os.pathsep):
        fpath = os.path.join(path, filename)
        if is_executable(fpath):
            return fpath

    return None


def executable_exists(filename):
    """Returns true if the filename refers to an executable file,
    either by relative or full path, or if the executable is found
    on the current PATH."""
    exec_path = which_executable(filename)

    return exec_path and is_executable(exec_path)


def missing_executables(filenames):
    result = []
    for filename in safe_coerce_to_frozenset(filenames):
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
    _sh_wrapper(shutil.move, source, destination)


def copy_file(source, destination):
    """Wrapper around shutils which ensures that the
    destination directory exists before copying the file."""
    _sh_wrapper(shutil.copy, source, destination)


def open_ro(filename):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(filename)
    try:
        header = handle.read(2)

        if header == "\x1f\x8b":
            handle.close()
            # TODO: Re-use handle (fileobj)
            handle = gzip.open(filename)
        elif header == "BZ":
            handle.close()
            handle = bz2.BZ2File(filename)
        else:
            handle.seek(0)

        return handle
    except:
        handle.close()
        raise


def try_remove(filename):
    """Tries to remove a file. Unlike os.remove, the function does not
    raise an exception if the file does not exist, but does raise
    exceptions on other errors. The return value reflects whether or
    not the file was actually removed."""
    return _try_rm_wrapper(os.remove, filename)


def try_rmdir(filename):
    """Tries to remove a directory. Unlike os.rmdir, the function does not raise
    an exception if the file does not exist, but does raise exceptions on other
    errors. The return value reflects whether or not the file was actually
    removed."""
    return _try_rm_wrapper(os.rmdir, filename)


def try_rmtree(filename):
    """Tries to remove a dir-tree. Unlike shutil.rmtree, the function does not raise
    an exception if the file does not exist, but does raise exceptions on other
    errors. The return value reflects whether or not the file was actually
    removed."""
    return _try_rm_wrapper(shutil.rmtree, filename)


def describe_files(files):
    """Return a text description of a set of files."""
    files = _validate_filenames(files)

    if not files:
        return "No files"
    elif len(files) == 1:
        return repr(files[0])

    glob_files = _get_files_glob(files, max_differences = 2)
    if glob_files:
        return repr(glob_files)

    paths = set(os.path.dirname(filename) for filename in files)
    if len(paths) == 1:
        return "%i files in '%s'" % (len(files), paths.pop())
    return "%i files" % (len(files),)


def describe_paired_files(files_1, files_2):
    """Return a text description of a set of paired filenames; the
    sets must be of the same length, a the description will depend
    on the overall similarity of the filenames / paths. If 'files_2'
    is empty, this function is the equivalent of calling
    'describe_files' with 'files_1' as the argument. In all other
    cases the length of the two sets must be the same."""
    files_1 = _validate_filenames(files_1)
    files_2 = _validate_filenames(files_2)

    if files_1 and not files_2:
        return describe_files(files_1)
    elif len(files_1) != len(files_2):
        raise ValueError("Unequal number of files for mate 1 vs mate 2 reads: %i vs %i" \
                         % (len(files_1), len(files_2)))

    glob_files_1 = _get_files_glob(files_1, 3)
    glob_files_2 = _get_files_glob(files_2, 3)
    if glob_files_1 and glob_files_2:
        final_glob = _get_files_glob((glob_files_1, glob_files_2), 1, show_differences = True)
        if final_glob:
            return repr(final_glob)

    fnames = files_1 + files_2
    paths = set(os.path.dirname(fname) for fname in fnames)
    if len(paths) == 1:
        return "%i pair(s) of files in '%s'" % (len(files_1), paths.pop())
    return "%i pair(s) of files" % (len(files_1),)


def _get_files_glob(filenames, max_differences = 1, show_differences = False):
    """Tries to generate a glob-string for a set of filenames, containing
    at most 'max_differences' different columns. If more differences are
    found, or if the length of filenames vary, None is returned."""
    # File lengths must be the same, otherwise we'd have to do MSA
    if len(set(map(len, filenames))) > 1:
        return None

    glob_fname, differences = [], 0
    for chars in zip(*filenames):
        if "?" in chars:
            chars = ('?',)

        if len(frozenset(chars)) > 1:
            if show_differences:
                chars = ("[%s]" % ("".join(sorted(chars))),)
            else:
                chars = ("?",)
            differences += 1
        glob_fname.append(chars[0])

    if differences > max_differences:
        return None

    return "".join(glob_fname)


def _validate_filenames(filenames):
    """Sanity checks for filenames handled by
    'describe_files' and 'describe_paired_files."""
    filenames = safe_coerce_to_tuple(filenames)
    for filename in filenames:
        if not isinstance(filename, types.StringTypes):
            raise ValueError("Only string types are allowed for filenames, not %s" \
                             % (filename.__class__.__name__,))
    return filenames


def _sh_wrapper(func, source, destination):
    """Runs an 'shutil' function ('func') which takes an 'source' and
    a 'destination' argument (e.g. copy/move/etc.), but silently
    handles the case where the destination directory does not exist.

    If this is the case, the function will first create the destination
    directory, and then retry the function."""
    try:
        func(source, destination)
    except IOError, error:
        if (error.errno == errno.ENOENT):
            if source and destination and os.path.exists(source):
                dirname = os.path.dirname(destination)
                make_dirs(dirname)
                func(source, destination)
                return
        elif (error.errno == errno.ENOSPC):
            # Not enough space; remove partial file
            os.unlink(destination)
        raise


def _try_rm_wrapper(func, fpath):
    """Takes a function (e.g. os.remove / os.rmdir), and attempts to remove a
    path; returns true if that path was succesfully remove, and false if it did
    not exist."""
    try:
        func(fpath)
        return True
    except OSError, error:
        if error.errno != errno.ENOENT:
            raise
        return False

