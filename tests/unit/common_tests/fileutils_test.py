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
import errno

import nose
from nose.tools import assert_equal, assert_raises

import pypeline
from tests.common.utils import with_temp_folder, monkeypatch, set_cwd, \
     set_file_contents, \
     get_file_contents, \
     assert_in

from pypeline.common.fileutils import \
     add_postfix, \
     swap_ext, \
     reroot_path, \
     create_temp_dir, \
     missing_files, \
     modified_after, \
     is_executable, \
     executable_exists, \
     missing_executables, \
     make_dirs, \
     move_file, \
     copy_file, \
     open_ro, \
     try_remove, \
     describe_files


################################################################################
################################################################################
## Tests for 'add_postfix'

def test_add_postfix__no_postfix():
    assert_equal(add_postfix("name.foo", ""), "name.foo")

def test_add_postfix__dot_postfix():
    assert_equal(add_postfix("name.foo", ".pf"), "name.pf.foo")

def test_add_postfix__underscore_postfix():
    assert_equal(add_postfix("name.foo", "_pf"), "name_pf.foo")


def test_add_postfix__no_ext__no_postfix():
    assert_equal(add_postfix("name", ""), "name")

def test_add_postfix__no_ext__dot_postfix():
    assert_equal(add_postfix("name", ".pf"), "name.pf")

def test_add_postfix__no_ext__underscore_postfix():
    assert_equal(add_postfix("name", "_pf"), "name_pf")




################################################################################
################################################################################
## Tests for 'swap_ext'

def test_swap_ext__has_ext_vs_empty_ext():
    assert_equal(swap_ext("name.foo", ""), "name")

def test_swap_ext__empty_ext_vs_empty_ext():
    assert_equal(swap_ext("name", ""), "name")

def test_swap_ext__has_ext_vs_dot_ext():
    assert_equal(swap_ext("name.foo", "."), "name")

def test_swap_ext__dot_ext_vs_dot_ext():
    assert_equal(swap_ext("name.", "."), "name")


def test_swap_ext__multiple__has_ext_vs_empty_ext():
    assert_equal(swap_ext("name.foo.bar", ""), "name.foo")

def test_swap_ext__multiple__has_ext_vs_dot_ext():
    assert_equal(swap_ext("name.foo.bar", "."), "name.foo")

def test_swap_ext__multiple__dot_ext_vs_dot_ext():
    assert_equal(swap_ext("name.foo.", "."), "name.foo")


def test_swap_ext__has_ext_vs_new_ext():
    assert_equal(swap_ext("name.foo", "bar"), "name.bar")

def test_swap_ext__has_ext_vs_new_dot_ext():
    assert_equal(swap_ext("name.foo", ".bar"), "name.bar")

def test_swap_ext__empty_ext_vs_new_ext():
    assert_equal(swap_ext("name", "bar"), "name.bar")

def test_swap_ext__dot_ext_vs_new_dot_ext():
    assert_equal(swap_ext("name", ".bar"), "name.bar")



################################################################################
################################################################################
## Tests for 'reroot_path'

def test_reroot_path__empty_root():
    assert_equal(reroot_path("", "/etc/apt/sources.list"), "sources.list")

def test_reroot_path__empty_path():
    assert_equal(reroot_path("/etc/apt", ""), "/etc/apt/")


def test_reroot_path__abs_abs__wo_final_dash():
    assert_equal(reroot_path("/etc/apt", "/tmp/sources.list"), "/etc/apt/sources.list")

def test_reroot_path__abs_abs__w_final_dash():
    assert_equal(reroot_path("/etc/apt/", "/tmp/sources.list"), "/etc/apt/sources.list")


def test_reroot_path__abs_rel__wo_final_dash():
    assert_equal(reroot_path("/etc/apt", "tmp/sources.list"), "/etc/apt/sources.list")

def test_reroot_path__abs_rel__w_final_dash():
    assert_equal(reroot_path("/etc/apt/", "tmp/sources.list"), "/etc/apt/sources.list")


def test_reroot_path__rel_abs__wo_final_dash():
    assert_equal(reroot_path("etc/apt", "/tmp/sources.list"), "etc/apt/sources.list")

def test_reroot_path__rel_abs__w_final_dash():
    assert_equal(reroot_path("etc/apt/", "/tmp/sources.list"), "etc/apt/sources.list")


def test_reroot_path__rel_rel__wo_final_dash():
    assert_equal(reroot_path("etc/apt", "tmp/sources.list"), "etc/apt/sources.list")

def test_reroot_path__rel_rel__w_final_dash():
    assert_equal(reroot_path("etc/apt/", "tmp/sources.list"), "etc/apt/sources.list")




################################################################################
################################################################################
## Tests for 'create_temp_dir'

@with_temp_folder
def test_create_temp_dir__create(temp_folder):
    tmp_dir_1 = create_temp_dir(temp_folder)
    tmp_dir_2 = create_temp_dir(temp_folder)
    assert os.path.exists(tmp_dir_1)
    assert os.path.exists(tmp_dir_2)

@with_temp_folder
def test_create_temp_dir__empty(temp_folder):
    tmp_dir  = create_temp_dir(temp_folder)
    contents = os.listdir(tmp_dir)
    assert not contents

@with_temp_folder
def test_create_temp_dir__permissions(temp_folder):
    tmp_dir = create_temp_dir(temp_folder)
    stats   = os.stat(tmp_dir)
    assert_equal(stats.st_mode & 0777, 0700)

@with_temp_folder
def test_create_temp_dir__creation_preempted(temp_folder):
    unwrapped, preempted_once = os.makedirs, []
    def _wrap_os_makedirs(*args, **kwargs):
        # Simulate somebody else creating the directory first
        if not preempted_once:
            unwrapped(*args, **kwargs)
            preempted_once.append(True)
        unwrapped(*args, **kwargs)

    with monkeypatch("os.makedirs", _wrap_os_makedirs):
        assert not os.listdir(temp_folder)
        work_dir = create_temp_dir(temp_folder)
        assert os.path.exists(temp_folder)
        dirs = os.listdir(temp_folder)
        assert_equal(len(dirs), 2)
        assert_in(os.path.basename(work_dir), dirs)
        assert bool(preempted_once)

def test_create_temp_dir__permission_denied():
    def _wrap_os_makedirs(*_args, **_kwargs):
        raise OSError((errno.EACCES, "Simulated premission denied"))

    with monkeypatch("os.makedirs", _wrap_os_makedirs):
        assert_raises(OSError, create_temp_dir, "/tmp")


################################################################################
################################################################################
## Tests for 'missing_files'

def test_missing_files__file_exists():
    assert_equal(missing_files(["tests/data/empty_file_1"]), [])

def test_missing_files__file_doesnt_exist():
    assert_equal(missing_files(["tests/data/missing_file_1"]),
                 ["tests/data/missing_file_1"])

def test_missing_files__mixed_files():
    files = ["tests/data/missing_file_1",
             "tests/data/empty_file_1"]
    result = ["tests/data/missing_file_1"]

    assert_equal(missing_files(files), result)




################################################################################
################################################################################
## Tests for 'modified_after'

def test_modified_after__modified_after():
    assert modified_after("tests/data/timestamp_a_younger", "tests/data/timestamp_a_older")
    assert modified_after("tests/data/timestamp_a_younger", "tests/data/timestamp_b_older")
    assert modified_after("tests/data/timestamp_b_younger", "tests/data/timestamp_a_older")

def test_modified_after__not_modified_after():
    assert not modified_after("tests/data/timestamp_a_older", "tests/data/timestamp_a_younger")
    assert not modified_after("tests/data/timestamp_a_older", "tests/data/timestamp_b_younger")
    assert not modified_after("tests/data/timestamp_b_older", "tests/data/timestamp_a_younger")

def test_modified_after__same_file():
    assert not modified_after("tests/data/timestamp_a_older", "tests/data/timestamp_a_older")
    assert not modified_after("tests/data/timestamp_a_older", "tests/data/timestamp_b_older")
    assert not modified_after("tests/data/timestamp_b_older", "tests/data/timestamp_a_older")




################################################################################
################################################################################
## Tests for 'is_executable'

def test_is_executable__full_path__is_executable():
    assert is_executable("/bin/ls")

def test_is_executable__full_path__is_non_executable():
    assert not is_executable("/etc/fstab")

def test_is_executable__rel_path__is_executable():
    assert is_executable("tests/unit/run")

def test_is_executable__rel_path__is_non_executable():
    assert not is_executable("tests/data/empty_file_1")




################################################################################
################################################################################
## Tests for 'executable_exists'

def test_executable_exists__executable():
    assert executable_exists("ls")

def test_executable_exists__non_executable():
    assert not executable_exists("lsxxxx")

def test_executable_exists__full_path__is_executable():
    assert executable_exists("/bin/ls")

def test_executable_exists__full_path__is_non_executable():
    assert not executable_exists("/etc/fstab")

def test_executable_exists__rel_path__is_executable():
    assert executable_exists("tests/unit/run")

def test_executable_exists__rel_path__is_non_executable():
    assert not executable_exists("tests/data/empty_file_1")




################################################################################
################################################################################
## Tests for 'missing_executables'

def test_missing_executables__executable():
    assert_equal(missing_executables(["ls"]), [])

def test_missing_executables__non_executable():
    assert_equal(missing_executables(["lsxxxx"]), ["lsxxxx"])

def test_missing_executables__mixed():
    assert_equal(missing_executables(["lsxxxx", "ls"]), ["lsxxxx"])




################################################################################
################################################################################
## Tests for 'make_dirs'

@with_temp_folder
def test_make_dirs__create_dir(temp_folder):
    assert not os.listdir(temp_folder)
    assert make_dirs(os.path.join(temp_folder, "test123"))
    assert_equal(os.listdir(temp_folder), ["test123"])

@with_temp_folder
def test_make_dirs__return_values(temp_folder):
    assert make_dirs(os.path.join(temp_folder, "test234"))
    assert not make_dirs(os.path.join(temp_folder, "test234"))

@with_temp_folder
def test_make_dirs__subdirs_return_values(temp_folder):
    assert make_dirs(os.path.join(temp_folder, "test"))
    assert make_dirs(os.path.join(temp_folder, "test", "234"))
    assert not make_dirs(os.path.join(temp_folder, "test", "234"))

@with_temp_folder
def test_make_dirs__sub_directories(temp_folder):
    assert not os.listdir(temp_folder)
    assert make_dirs(os.path.join(temp_folder, "test", "123"))
    assert_equal(os.listdir(temp_folder), ["test"])
    assert_equal(os.listdir(os.path.join(temp_folder, "test")), ["123"])

@with_temp_folder
def test_make_dirs__permissions(temp_folder):
    work_dir = os.path.join(temp_folder, "test_1")
    assert make_dirs(work_dir, mode = 0511)
    stats   = os.stat(work_dir)
    assert_equal(oct(stats.st_mode & 0777), oct(0511))

@with_temp_folder
def test_make_dirs__creation_preemted(temp_folder):
    unwrapped, preempted = os.makedirs, []
    def _wrap_os_makedirs(*args, **kwargs):
        # Simulate somebody else creating the directory first
        preempted.append(True)
        unwrapped(*args, **kwargs)
        unwrapped(*args, **kwargs)

    with monkeypatch("os.makedirs", _wrap_os_makedirs):
        work_folder = os.path.join(temp_folder, "test")
        assert not make_dirs(work_folder)
        assert os.path.exists(work_folder)
        assert_equal(os.listdir(temp_folder), ["test"])
        assert_equal(preempted, [True])

@nose.tools.raises(ValueError)
def test_make_dirs__empty_directory():
    make_dirs("")




################################################################################
################################################################################
## Tests for 'move_file'

@with_temp_folder
def test_move_file__simple_move(temp_folder):
    file_1 = os.path.join(temp_folder, "file_1")
    file_2 = os.path.join(temp_folder, "file_2")
    assert_equal(os.listdir(temp_folder), [])
    set_file_contents(file_1, "1")
    assert_equal(os.listdir(temp_folder), ["file_1"])
    move_file(file_1, file_2)
    assert_equal(os.listdir(temp_folder), ["file_2"])
    assert_equal(get_file_contents(file_2), "1")

@with_temp_folder
def test_move_file__simple_move_in_cwd(temp_folder):
    with set_cwd(temp_folder):
        assert_equal(os.listdir("."), [])
        set_file_contents("file_1", "1")
        assert_equal(os.listdir("."), ["file_1"])
        move_file("file_1", "file_2")
        assert_equal(os.listdir("."), ["file_2"])
        assert_equal(get_file_contents("file_2"), "1")


@with_temp_folder
def test_move_file__move_to_existing_folder(temp_folder):
    assert make_dirs(os.path.join(temp_folder, "src"))
    assert make_dirs(os.path.join(temp_folder, "dst"))
    file_1 = os.path.join(temp_folder, "src", "file_1")
    file_2 = os.path.join(temp_folder, "dst", "file_2")
    set_file_contents(file_1, "2")
    move_file(file_1, file_2)
    assert_equal(os.listdir(os.path.dirname(file_1)), [])
    assert_equal(os.listdir(os.path.dirname(file_2)), ["file_2"])
    assert_equal(get_file_contents(file_2), "2")

@with_temp_folder
def test_move_file__move_to_new_folder(temp_folder):
    assert make_dirs(os.path.join(temp_folder, "src"))
    file_1 = os.path.join(temp_folder, "src", "file_1")
    file_2 = os.path.join(temp_folder, "dst", "file_2")
    set_file_contents(file_1, "2")
    move_file(file_1, file_2)
    assert_equal(os.listdir(os.path.dirname(file_1)), [])
    assert_equal(os.listdir(os.path.dirname(file_2)), ["file_2"])
    assert_equal(get_file_contents(file_2), "2")

@with_temp_folder
def test_move_file__move_to_different_folder(temp_folder):
    with set_cwd(temp_folder):
        set_file_contents("file_1", "3")
        move_file("file_1", "dst/file_1")
        assert_equal(os.listdir("."), ["dst"])
        assert_equal(os.listdir("dst"), ["file_1"])
        assert_equal(get_file_contents("dst/file_1"), "3")

@with_temp_folder
def test_move_file__overwrite(temp_folder):
    with set_cwd(temp_folder):
        set_file_contents("file_1", "4")
        set_file_contents("file_2", "5")
        move_file("file_1", "file_2")
        assert_equal(os.listdir("."), ["file_2"])
        assert_equal(get_file_contents("file_2"), "4")




################################################################################
################################################################################
## Tests for 'copy_file'

@with_temp_folder
def test_copy_file__simple_copy(temp_folder):
    file_1 = os.path.join(temp_folder, "file_1")
    file_2 = os.path.join(temp_folder, "file_2")
    assert_equal(os.listdir(temp_folder), [])
    set_file_contents(file_1, "1")
    assert_equal(os.listdir(temp_folder), ["file_1"])
    copy_file(file_1, file_2)
    assert_equal(set(os.listdir(temp_folder)), set(["file_1", "file_2"]))
    assert_equal(get_file_contents(file_1), "1")
    assert_equal(get_file_contents(file_2), "1")

@with_temp_folder
def test_copy_file__simple_copy_in_cwd(temp_folder):
    with set_cwd(temp_folder):
        assert_equal(os.listdir("."), [])
        set_file_contents("file_1", "1")
        assert_equal(os.listdir("."), ["file_1"])
        copy_file("file_1", "file_2")
        assert_equal(set(os.listdir(".")), set(["file_1", "file_2"]))
        assert_equal(get_file_contents("file_1"), "1")
        assert_equal(get_file_contents("file_2"), "1")

@with_temp_folder
def test_copy_file__copy_to_existing_folder(temp_folder):
    assert make_dirs(os.path.join(temp_folder, "src"))
    assert make_dirs(os.path.join(temp_folder, "dst"))
    file_1 = os.path.join(temp_folder, "src", "file_1")
    file_2 = os.path.join(temp_folder, "dst", "file_2")
    set_file_contents(file_1, "2")
    copy_file(file_1, file_2)
    assert_equal(os.listdir(os.path.dirname(file_1)), ["file_1"])
    assert_equal(os.listdir(os.path.dirname(file_2)), ["file_2"])
    assert_equal(get_file_contents(file_1), "2")
    assert_equal(get_file_contents(file_2), "2")

@with_temp_folder
def test_copy_file__copy_to_new_folder(temp_folder):
    assert make_dirs(os.path.join(temp_folder, "src"))
    file_1 = os.path.join(temp_folder, "src", "file_1")
    file_2 = os.path.join(temp_folder, "dst", "file_2")
    set_file_contents(file_1, "2")
    copy_file(file_1, file_2)
    assert_equal(os.listdir(os.path.dirname(file_1)), ["file_1"])
    assert_equal(os.listdir(os.path.dirname(file_2)), ["file_2"])
    assert_equal(get_file_contents(file_1), "2")
    assert_equal(get_file_contents(file_2), "2")

@with_temp_folder
def test_copy_file__copy_to_different_folder(temp_folder):
    with set_cwd(temp_folder):
        set_file_contents("file_1", "3")
        copy_file("file_1", "dst/file_1")
        assert_equal(set(os.listdir(".")), set(["file_1", "dst"]))
        assert_equal(os.listdir("dst"), ["file_1"])
        assert_equal(get_file_contents("file_1"), "3")
        assert_equal(get_file_contents("dst/file_1"), "3")

@with_temp_folder
def test_copy_file__overwrite(temp_folder):
    with set_cwd(temp_folder):
        set_file_contents("file_1", "4")
        set_file_contents("file_2", "5")
        copy_file("file_1", "file_2")
        assert_equal(set(os.listdir(".")), set(["file_1", "file_2"]))
        assert_equal(get_file_contents("file_1"), "4")
        assert_equal(get_file_contents("file_2"), "4")


################################################################################
################################################################################
## Tests for 'open'


def test_open_ro__uncompressed():
    handle = open_ro('tests/data/fasta_file.fasta')
    try:
        assert_equal(handle.read(), b'>This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n')
    finally:
        handle.close()

def test_open_ro__gz():
    handle = open_ro('tests/data/fasta_file.fasta.gz')
    try:
        assert_equal(handle.read(), b'>This_is_GZipped_FASTA!\nACGTN\n>This_is_ALSO_GZipped_FASTA!\nCGTNA\n')
    finally:
        handle.close()

def test_open_ro__bz2():
    handle = open_ro('tests/data/fasta_file.fasta.bz2')
    try:
        assert_equal(handle.read(), b'>This_is_BZ_FASTA!\nCGTNA\n>This_is_ALSO_BZ_FASTA!\nACGTN\n')
    finally:
        handle.close()


class OddException(RuntimeError):
    pass

def test_open_ro__close_handle_on_error():
    class _FileMock:
        def __init__(self, filename):
            self._close_called = False
            assert_equal(filename, "/var/abc")
        def read(self, *_args, **_kwargs):
            # pylint: disable=R0201
            raise OddException("ARGH!")
        def close(self):
            self._close_called = True

    try:
        pypeline.common.fileutils.open = _FileMock
        assert_raises(OddException, open_ro, "/var/abc")
    finally:
        del pypeline.common.fileutils.open



################################################################################
################################################################################
## Tests for 'try_remove'

@with_temp_folder
def test_try_remove(temp_folder):
    fpath = os.path.join(temp_folder, "test.txt")
    set_file_contents(fpath, "1 2 3")
    assert try_remove(fpath)


@with_temp_folder
def test_try_remove__missing(temp_folder):
    fpath = os.path.join(temp_folder, "test.txt")
    assert not try_remove(fpath)


@with_temp_folder
@nose.tools.raises(OSError)
def test_try_remove__non_file(temp_folder):
    try_remove(temp_folder)



################################################################################
################################################################################
## Tests for 'describe_files'

def test_describe_files__no_files():
    assert_equal(describe_files(()), "No files")

def test_describe_files__single_file():
    fpath = "/var/foo/bar"
    assert_equal(describe_files((fpath,)), repr(fpath))

def test_describe_files__same_path_abs():
    fpaths = ("/var/foo/bar", "/var/foo/foo")
    assert_equal(describe_files(fpaths), "2 files in '/var/foo'")

def test_describe_files__different_paths_abs():
    fpaths = ("/var/foo/bar", "/var/bar/foo")
    assert_equal(describe_files(fpaths), "2 files")

def test_describe_files__same_path_rel():
    fpaths = ("var/foo/bar", "var/foo/foo")
    assert_equal(describe_files(fpaths), "2 files in 'var/foo'")

def test_describe_files__different_paths_rel():
    fpaths = ("var/foo/bar", "var/bar/foo")
    assert_equal(describe_files(fpaths), "2 files")

def test_describe_files__iterable():
    fpaths = iter(("/var/foo/bar", "/var/foo/foo"))
    assert_equal(describe_files(fpaths), "2 files in '/var/foo'")
