#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import bz2
import errno
import gzip
import io
import os
import random
import shutil
import tempfile
from datetime import datetime
from os import fspath
from typing import (
    IO,
    Any,
    Callable,
    Iterable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

from .utilities import safe_coerce_to_tuple

PathTypes = Union[str, "os.PathLike[str]"]


def add_postfix(filename: PathTypes, postfix: str) -> str:
    """Ads a postfix to a filename (before any extensions that filename may have)."""
    filename, ext = os.path.splitext(filename)
    return filename + postfix + ext


def swap_ext(filename: str, ext: str) -> str:
    """Replaces the existing extension of a filename with the specified extension,
    ensuring that the extension is prefixed by a '.'. If no extension is specified,
    other than potentially a dot, the existing extension is stripped."""
    filename, _ = os.path.splitext(filename)
    if ext in ("", "."):
        return filename

    if not ext.startswith("."):
        ext = "." + ext

    return filename + ext


def reroot_path(root: PathTypes, filename: str) -> str:
    """Returns the basename of filename, joined to root."""
    return os.path.join(root, os.path.basename(filename))


def create_temp_dir(root: PathTypes) -> str:
    """Creates a temporary directory, accessible only by the owner at the specified
    location. The folder name is includes the current time and includes a random
    component. Only the current user has access.
    """
    now = datetime.now()  # noqa: DTZ005

    return tempfile.mkdtemp(
        prefix=now.strftime("%Y%m%d_%H%M%S_"),
        dir=root,
    )


    for _ in range(10_000):
        path = _generate_path()
        if make_dirs(path, mode=0o700):
            return path

    raise FileExistsError(errno.EEXIST, "Could not create new temp directory")


def missing_files(filenames: Iterable[PathTypes]) -> List[str]:
    """Given a list of filenames, returns a list of those that
    does not exist. Note that this function does not differentiate
    between files and folders."""
    missing: List[str] = []
    for filename in safe_coerce_to_tuple(filenames):
        filename = fspath(filename)
        if not os.path.exists(filename):
            missing.append(filename)
    return missing


def missing_executables(filenames: Iterable[PathTypes]) -> list[str]:
    missing: list[str] = []
    for filename in filenames:
        value = fspath(filename)
        if not shutil.which(value):
            missing.append(value)
    return missing


def make_dirs(directory: PathTypes, mode: int = 0o777) -> bool:
    """Wrapper around os.makedirs to make it suitable for using
    in a multithreaded/multiprocessing environment: Unlike the
    regular function, this wrapper does not throw an exception if
    the directory already exists, which may happen if another
    thread/process created the directory during the function call.

    Returns true if a new directory was created, false if it
    already existed. Other errors result in exceptions."""
    if not directory:
        raise ValueError("Empty directory passed to make_dirs()")

    try:
        os.makedirs(fspath(directory), mode=mode)
        return True
    except OSError as error:
        # make_dirs be called by multiple subprocesses at the same time,
        # so only raise if the actual creation of the folder failed
        if error.errno != errno.EEXIST:
            raise
        return False


def move_file(source: PathTypes, destination: PathTypes) -> None:
    """Wrapper around shutils which ensures that the
    destination directory exists before moving the file.

    In addition, the (partial) destination file is automatically removed if the copy
    fails with an out-of-space error.
    """
    _sh_wrapper(_atomic_file_move, source, destination)


def copy_file(source: PathTypes, destination: PathTypes) -> None:
    """Wrapper around shutils which ensures that the
    destination directory exists before copying the file.

    In addition, the (partial) destination file is automatically removed if the copy
    fails with an out-of-space error.
    """
    _sh_wrapper(_atomic_file_copy, source, destination)


def open_rb(filename: PathTypes) -> IO[bytes]:
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(fspath(filename), "rb")

    try:
        header = handle.peek(2)

        if header.startswith(b"\x1f\x8b"):
            return cast(IO[bytes], _GzipFile(mode="rb", fileobj=handle))
        elif header.startswith(b"BZ"):
            return _BZ2File(handle, "rb")
        else:
            return handle
    except:
        handle.close()
        raise


def open_rt(filename: PathTypes) -> IO[str]:
    return io.TextIOWrapper(open_rb(filename))


def try_remove(filename: PathTypes) -> bool:
    """Tries to remove a file. Unlike os.remove, the function does not
    raise an exception if the file does not exist, but does raise
    exceptions on other errors. The return value reflects whether or
    not the file was actually removed."""
    return _try_rm_wrapper(os.remove, filename)


def try_rmtree(filename: PathTypes) -> bool:
    """Tries to remove a dir-tree. Unlike shutil.rmtree, the function does not raise
    an exception if the file does not exist, but does raise exceptions on other
    errors. The return value reflects whether or not the file was actually
    removed."""
    return _try_rm_wrapper(shutil.rmtree, filename)


def try_rmdirs(dirpath: PathTypes) -> bool:
    try:
        items = os.scandir(dirpath)
    except FileNotFoundError:
        return True

    is_empty = True
    for it in items:
        if it.is_symlink() or not it.is_dir():
            is_empty = False
        elif not try_rmdirs(it.path):
            is_empty = False

    if is_empty:
        # If _try_rm_wrapper returns then the folder no longer exists
        _try_rm_wrapper(os.rmdir, dirpath)

    return is_empty


def describe_files(files: Iterable[str]) -> str:
    """Return a text description of a set of files."""
    files = validate_filenames(files)

    if not files:
        return "No files"
    elif len(files) == 1:
        return repr(files[0])

    glob_files = get_files_glob(files, max_differences=2)
    if glob_files:
        return repr(glob_files)

    paths = set(os.path.dirname(filename) for filename in files)
    if len(paths) == 1:
        return "%i files in '%s'" % (len(files), paths.pop())
    return "%i files" % (len(files),)


def describe_paired_files(files_1: Iterable[str], files_2: Iterable[str]) -> str:
    """Return a text description of a set of paired filenames; the
    sets must be of the same length, a the description will depend
    on the overall similarity of the filenames / paths. If 'files_2'
    is empty, this function is the equivalent of calling
    'describe_files' with 'files_1' as the argument. In all other
    cases the length of the two sets must be the same."""
    files_1 = validate_filenames(files_1)
    files_2 = validate_filenames(files_2)

    if files_1 and not files_2:
        return describe_files(files_1)
    elif len(files_1) != len(files_2):
        raise ValueError(
            "Unequal number of files for mate 1 vs mate 2 reads: %i vs %i"
            % (len(files_1), len(files_2))
        )

    glob_files_1 = get_files_glob(files_1, 3)
    glob_files_2 = get_files_glob(files_2, 3)
    if glob_files_1 and glob_files_2:
        final_glob = get_files_glob(
            (glob_files_1, glob_files_2), 1, show_differences=True
        )
        if final_glob:
            return repr(final_glob)

    paths = {os.path.dirname(fname) for fname in (files_1 + files_2)}
    if len(paths) == 1:
        return "%i pair(s) of files in '%s'" % (len(files_1), paths.pop())
    return "%i pair(s) of files" % (len(files_1),)


def get_files_glob(
    filenames: Sequence[str],
    max_differences: int = 1,
    show_differences: bool = False,
) -> Optional[str]:
    """Tries to generate a glob-string for a set of filenames, containing
    at most 'max_differences' different columns. If more differences are
    found, or if the length of filenames vary, None is returned."""
    # File lengths must be the same, otherwise we'd have to do MSA
    if len(set(map(len, filenames))) > 1:
        return None

    glob_fname: List[str] = []
    differences = 0
    for chars in zip(*filenames):
        if "?" in chars:
            chars = ("?",)

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


def validate_filenames(filenames: Iterable[str]) -> Tuple[str, ...]:
    """Sanity checks for filenames handled by
    'describe_files' and 'describe_paired_files."""
    return tuple(fspath(filename) for filename in safe_coerce_to_tuple(filenames))


def _atomic_file_move(source: PathTypes, destination: PathTypes):
    try:
        return os.rename(source, destination)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise

    # Move across filesystem/device
    if os.path.islink(source):
        os.symlink(os.readlink(source), destination)
        os.unlink(source)
        return

    # Copy with copy2 to preserve metadata
    _atomic_file_copy(source, destination, shutil.copy2)
    os.unlink(source)


def _atomic_file_copy(
    source: PathTypes,
    destination: PathTypes,
    func: Optional[Callable[[PathTypes, PathTypes], Any]] = None,
):
    # Ensure that hard failures during copy does not leave anything at destination
    postfix = "{:04x}".format(random.getrandbits(16))
    temp_destination = "{}.{}.tmp".format(destination, postfix)
    # Allow function to be monkey-patched for testing
    func = shutil.copy if func is None else func

    try:
        func(source, temp_destination)
    except OSError as error:
        if error.errno == errno.ENOSPC:
            # Not enough space; remove partial file
            os.unlink(temp_destination)
        raise

    os.rename(temp_destination, destination)


def _sh_wrapper(
    func: Callable[[PathTypes, PathTypes], Any],
    source: PathTypes,
    destination: PathTypes,
) -> None:
    """Runs an 'shutil' function ('func') which takes an 'source' and
    a 'destination' argument (e.g. copy/move/etc.), but silently
    handles the case where the destination directory does not exist.

    If this is the case, the function will first create the destination
    directory, and then retry the function."""
    source = fspath(source)
    destination = fspath(destination)

    try:
        func(source, destination)
    except FileNotFoundError:
        if source and destination and os.path.exists(source):
            make_dirs(os.path.dirname(destination))
            return func(source, destination)
        raise


def _try_rm_wrapper(func: Callable[[Any], Any], fpath: PathTypes) -> bool:
    """Takes a function (e.g. os.remove / os.rmdir), and attempts to remove a
    path; returns true if that path was successfully remove, and false if it did
    not exist."""
    try:
        func(fspath(fpath))
        return True
    except FileNotFoundError:
        return False


class _GzipFile(gzip.GzipFile):
    "Wrapper ensuring that passed filehandles are properly closed"

    def close(self) -> None:
        fileobj: Any = self.fileobj
        super().close()
        if hasattr(fileobj, "close"):
            fileobj.close()


class _BZ2File(bz2.BZ2File):
    "Wrapper ensuring that passed filehandles are properly closed"

    def close(self) -> None:
        fileobj = self._fp  # type: ignore
        super().close()
        if hasattr(fileobj, "close"):
            fileobj.close()
