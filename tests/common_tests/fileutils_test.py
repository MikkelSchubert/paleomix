#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import errno
import os
import stat

from pathlib import Path
from typing import Any
from unittest.mock import ANY, call, DEFAULT, Mock, patch

import pytest

from paleomix.common.testing import SetWorkingDirectory

from paleomix.common.fileutils import (
    add_postfix,
    swap_ext,
    reroot_path,
    create_temp_dir,
    missing_files,
    missing_executables,
    make_dirs,
    move_file,
    copy_file,
    open_ro,
    try_remove,
    try_rmtree,
    describe_files,
    describe_paired_files,
)


###############################################################################
###############################################################################
# Setup timestamps for test files


def test_dir() -> str:
    return os.path.dirname(os.path.dirname(__file__))


def test_file(*args: str) -> str:
    return os.path.join(test_dir(), "data", *args)


###############################################################################
###############################################################################
# Tests for 'add_postfix'


def test_add_postfix__no_postfix() -> None:
    assert add_postfix("name.foo", "") == "name.foo"


def test_add_postfix__dot_postfix() -> None:
    assert add_postfix("name.foo", ".pf") == "name.pf.foo"


def test_add_postfix__underscore_postfix() -> None:
    assert add_postfix("name.foo", "_pf") == "name_pf.foo"


def test_add_postfix__no_ext__no_postfix() -> None:
    assert add_postfix("name", "") == "name"


def test_add_postfix__no_ext__dot_postfix() -> None:
    assert add_postfix("name", ".pf") == "name.pf"


def test_add_postfix__no_ext__underscore_postfix() -> None:
    assert add_postfix("name", "_pf") == "name_pf"


###############################################################################
###############################################################################
# Tests for 'swap_ext'


def test_swap_ext__has_ext_vs_empty_ext() -> None:
    assert swap_ext("name.foo", "") == "name"


def test_swap_ext__empty_ext_vs_empty_ext() -> None:
    assert swap_ext("name", "") == "name"


def test_swap_ext__has_ext_vs_dot_ext() -> None:
    assert swap_ext("name.foo", ".") == "name"


def test_swap_ext__dot_ext_vs_dot_ext() -> None:
    assert swap_ext("name.", ".") == "name"


def test_swap_ext__multiple__has_ext_vs_empty_ext() -> None:
    assert swap_ext("name.foo.bar", "") == "name.foo"


def test_swap_ext__multiple__has_ext_vs_dot_ext() -> None:
    assert swap_ext("name.foo.bar", ".") == "name.foo"


def test_swap_ext__multiple__dot_ext_vs_dot_ext() -> None:
    assert swap_ext("name.foo.", ".") == "name.foo"


def test_swap_ext__has_ext_vs_new_ext() -> None:
    assert swap_ext("name.foo", "bar") == "name.bar"


def test_swap_ext__has_ext_vs_new_dot_ext() -> None:
    assert swap_ext("name.foo", ".bar") == "name.bar"


def test_swap_ext__empty_ext_vs_new_ext() -> None:
    assert swap_ext("name", "bar") == "name.bar"


def test_swap_ext__dot_ext_vs_new_dot_ext() -> None:
    assert swap_ext("name", ".bar") == "name.bar"


###############################################################################
###############################################################################
# Tests for 'reroot_path'


def test_reroot_path__empty_root() -> None:
    assert reroot_path("", "/etc/apt/sources.list") == "sources.list"


def test_reroot_path__empty_path() -> None:
    assert reroot_path("/etc/apt", "") == "/etc/apt/"


def test_reroot_path__abs_abs__wo_final_dash() -> None:
    assert reroot_path("/etc/apt", "/tmp/sources.list") == "/etc/apt/sources.list"


def test_reroot_path__abs_abs__w_final_dash() -> None:
    assert reroot_path("/etc/apt/", "/tmp/sources.list") == "/etc/apt/sources.list"


def test_reroot_path__abs_rel__wo_final_dash() -> None:
    assert reroot_path("/etc/apt", "tmp/sources.list") == "/etc/apt/sources.list"


def test_reroot_path__abs_rel__w_final_dash() -> None:
    assert reroot_path("/etc/apt/", "tmp/sources.list") == "/etc/apt/sources.list"


def test_reroot_path__rel_abs__wo_final_dash() -> None:
    assert reroot_path("etc/apt", "/tmp/sources.list") == "etc/apt/sources.list"


def test_reroot_path__rel_abs__w_final_dash() -> None:
    assert reroot_path("etc/apt/", "/tmp/sources.list") == "etc/apt/sources.list"


def test_reroot_path__rel_rel__wo_final_dash() -> None:
    assert reroot_path("etc/apt", "tmp/sources.list") == "etc/apt/sources.list"


def test_reroot_path__rel_rel__w_final_dash() -> None:
    assert reroot_path("etc/apt/", "tmp/sources.list") == "etc/apt/sources.list"


###############################################################################
###############################################################################
# Tests for 'create_temp_dir'


def test_create_temp_dir__create(tmp_path: Path) -> None:
    tmp_dir_1 = create_temp_dir(tmp_path)
    tmp_dir_2 = create_temp_dir(tmp_path)
    assert os.path.exists(tmp_dir_1)
    assert os.path.exists(tmp_dir_2)


def test_create_temp_dir__empty(tmp_path: Path) -> None:
    tmp_dir = create_temp_dir(tmp_path)
    contents = os.listdir(tmp_dir)
    assert not contents


def test_create_temp_dir__permissions(tmp_path: Path) -> None:
    tmp_dir = create_temp_dir(tmp_path)
    stats = os.stat(tmp_dir)
    assert stats.st_mode & 0 == 0


def test_create_temp_dir__creation_preempted(tmp_path: Path) -> None:
    makedirs = os.makedirs
    mock = Mock(wraps=makedirs)
    # Raise expection first, then pass second call to wrapped function
    mock.side_effect = [OSError(errno.EEXIST, "dir exists"), DEFAULT]
    with patch("os.makedirs", mock):
        work_dir = create_temp_dir(tmp_path)

    assert work_dir.startswith(str(tmp_path))

    # Must be called with different directories
    assert mock.mock_calls == [call(ANY, mode=0o750), call(ANY, mode=0o750)]
    call_1, call_2 = mock.mock_calls
    assert call_1[1] != call_2[1]


def test_create_temp_dir__permission_denied() -> None:
    with patch("os.makedirs") as mock:
        mock.side_effect = OSError(errno.EACCES, "Simulated premission denied")

        with pytest.raises(OSError, match="Simulated premission denied"):
            create_temp_dir("/tmp")


###############################################################################
###############################################################################
# Tests for 'missing_files'


def test_missing_files__file_exists() -> None:
    assert missing_files([test_file("empty_file_1")]) == []


def test_missing_files__file_doesnt_exist() -> None:
    assert missing_files([test_file("missing_file_1")]) == [test_file("missing_file_1")]


def test_missing_files__mixed_files() -> None:
    files = [test_file("missing_file_1"), test_file("empty_file_1")]
    result = [test_file("missing_file_1")]

    assert missing_files(files) == result


###############################################################################
###############################################################################
# Tests for 'missing_executables'


def test_missing_executables__executable() -> None:
    assert missing_executables(["ls"]) == []


def test_missing_executables__non_executable() -> None:
    assert missing_executables(["lsxxxx"]) == ["lsxxxx"]


def test_missing_executables__mixed() -> None:
    assert missing_executables(["lsxxxx", "ls"]) == ["lsxxxx"]


def test_missing_executables__with_path__missing(tmp_path: Path) -> None:
    executable = tmp_path / "ls"
    assert missing_executables([executable]) == [executable]


def test_missing_executables__with_path__exists(tmp_path: Path) -> None:
    executable = tmp_path / "ls"
    executable.touch()
    executable.chmod(0o700)

    assert missing_executables([executable]) == []


###############################################################################
###############################################################################
# Tests for 'make_dirs'


def test_make_dirs__create_dir(tmp_path: Path) -> None:
    assert not os.listdir(tmp_path)
    assert make_dirs(os.path.join(tmp_path, "test123"))
    assert os.listdir(tmp_path) == ["test123"]


def test_make_dirs__return_values(tmp_path: Path) -> None:
    assert make_dirs(os.path.join(tmp_path, "test234"))
    assert not make_dirs(os.path.join(tmp_path, "test234"))


def test_make_dirs__subdirs_return_values(tmp_path: Path) -> None:
    assert make_dirs(os.path.join(tmp_path, "test"))
    assert make_dirs(os.path.join(tmp_path, "test", "234"))
    assert not make_dirs(os.path.join(tmp_path, "test", "234"))


def test_make_dirs__sub_directories(tmp_path: Path) -> None:
    assert not os.listdir(tmp_path)
    assert make_dirs(os.path.join(tmp_path, "test", "123"))
    assert os.listdir(tmp_path) == ["test"]
    assert os.listdir(os.path.join(tmp_path, "test")) == ["123"]


def test_make_dirs__permissions(tmp_path: Path) -> None:
    work_dir = tmp_path / "test_1"
    assert make_dirs(work_dir, mode=0o511)
    stats = work_dir.stat()
    assert oct(stats.st_mode & 0o777) == oct(0o511)


def test_make_dirs__creation_preemted(tmp_path: Path) -> None:
    makedirs = os.makedirs

    def _wrap_os_makedirs(*args: Any, **kwargs: Any):
        # Simulate somebody else creating the directory first
        makedirs(*args, **kwargs)
        makedirs(*args, **kwargs)

    with patch("os.makedirs", _wrap_os_makedirs):
        work_folder = tmp_path / "test"
        assert not make_dirs(work_folder)
        assert os.path.exists(work_folder)
        assert os.listdir(tmp_path) == ["test"]


def test_make_dirs__permission_denied(tmp_path: Path) -> None:
    # Make temporary folder read-only
    mode = os.stat(tmp_path).st_mode
    ro_mode = mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    os.chmod(tmp_path, ro_mode)
    # Non OEXIST errors should be re-raised:
    with pytest.raises(OSError):
        make_dirs(os.path.join(tmp_path, "foo"))


def test_make_dirs__empty_directory() -> None:
    with pytest.raises(ValueError, match="Empty directory passed to make_dirs"):
        make_dirs("")


###############################################################################
###############################################################################
# Tests for 'move_file'


def test_move_file__simple_move(tmp_path: Path) -> None:
    file_1 = tmp_path / "file_1"
    file_2 = tmp_path / "file_2"
    assert os.listdir(tmp_path) == []
    file_1.write_text("1")
    assert os.listdir(tmp_path) == ["file_1"]
    move_file(file_1, file_2)
    assert os.listdir(tmp_path) == ["file_2"]
    assert file_2.read_text() == "1"


def test_move_file__simple_move_in_cwd(tmp_path: Path) -> None:
    with SetWorkingDirectory(tmp_path):
        assert os.listdir(".") == []
        (tmp_path / "file_1").write_text("1")
        assert os.listdir(".") == ["file_1"]
        move_file("file_1", "file_2")
        assert os.listdir(".") == ["file_2"]
        assert (tmp_path / "file_2").read_text() == "1"


def test_move_dirs__permission_denied(tmp_path: Path) -> None:
    dst_folder = tmp_path / "dst"
    file_1 = tmp_path / "file"
    file_2 = dst_folder / "file"
    file_1.write_text("1")

    # Make destination folder read-only
    assert make_dirs(tmp_path / "dst")
    mode = os.stat(dst_folder).st_mode
    ro_mode = mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    os.chmod(dst_folder, ro_mode)

    # Non ENOENT errors should be re-raised:
    with pytest.raises(IOError):
        move_file(file_1, file_2)


def test_move_file__move_to_existing_folder(tmp_path: Path) -> None:
    assert make_dirs(tmp_path / "src")
    assert make_dirs(tmp_path / "dst")
    file_1 = tmp_path / "src" / "file_1"
    file_2 = tmp_path / "dst" / "file_2"
    file_1.write_text("2")
    move_file(file_1, file_2)

    assert os.listdir(file_1.parent) == []
    assert os.listdir(file_2.parent) == ["file_2"]
    assert file_2.read_text() == "2"


def test_move_file__move_to_new_folder(tmp_path: Path) -> None:
    assert make_dirs(tmp_path / "src")
    file_1 = tmp_path / "src" / "file_1"
    file_2 = tmp_path / "dst" / "file_2"
    file_1.write_text("2")

    move_file(file_1, file_2)

    assert os.listdir(file_1.parent) == []
    assert os.listdir(file_2.parent) == ["file_2"]
    assert file_2.read_text() == "2"


def test_move_file__move_to_different_folder(tmp_path: Path) -> None:
    (tmp_path / "file_1").write_text("3")

    with SetWorkingDirectory(tmp_path):
        move_file("file_1", "dst/file_1")

    assert os.listdir(tmp_path), ["dst"]
    assert os.listdir(tmp_path / "dst"), ["file_1"]
    assert (tmp_path / "dst" / "file_1").read_text() == "3"


def test_move_file__overwrite(tmp_path: Path) -> None:
    (tmp_path / "file_1").write_text("4")
    (tmp_path / "file_2").write_text("5")

    with SetWorkingDirectory(tmp_path):
        move_file("file_1", "file_2")

    assert os.listdir(tmp_path) == ["file_2"]
    assert (tmp_path / "file_2").read_text() == "4"


def test_move_file__enoent_reraised_if_not_due_to_missing_folder() -> None:
    with pytest.raises(IOError):
        move_file("", "./dst")


###############################################################################
###############################################################################
# Tests for 'copy_file'


def test_copy_file__simple_copy(tmp_path: Path) -> None:
    file_1 = tmp_path / "file_1"
    file_2 = tmp_path / "file_2"
    file_1.write_text("1")
    assert os.listdir(tmp_path) == ["file_1"]
    copy_file(file_1, file_2)
    assert set(os.listdir(tmp_path)) == set(["file_1", "file_2"])
    assert file_1.read_text() == "1"
    assert file_2.read_text() == "1"


def test_copy_file__simple_copy_in_cwd(tmp_path: Path) -> None:
    file_1 = tmp_path / "file_1"
    file_2 = tmp_path / "file_2"
    file_1.write_text("1")

    with SetWorkingDirectory(tmp_path):
        copy_file("file_1", "file_2")

    assert set(os.listdir(tmp_path)) == set(["file_1", "file_2"])
    assert file_1.read_text() == "1"
    assert file_2.read_text() == "1"


def test_copy_file__copy_to_existing_folder(tmp_path: Path) -> None:
    assert make_dirs(tmp_path / "src")
    assert make_dirs(tmp_path / "dst")
    file_1 = tmp_path / "src" / "file_1"
    file_2 = tmp_path / "dst" / "file_2"
    file_1.write_text("2")
    copy_file(file_1, file_2)
    assert os.listdir(file_1.parent) == ["file_1"]
    assert os.listdir(file_2.parent) == ["file_2"]
    assert file_1.read_text() == "2"
    assert file_2.read_text() == "2"


def test_copy_file__copy_to_new_folder(tmp_path: Path) -> None:
    assert make_dirs(tmp_path / "src")
    file_1 = tmp_path / "src" / "file_1"
    file_2 = tmp_path / "dst" / "file_2"
    file_1.write_text("2")
    copy_file(file_1, file_2)
    assert os.listdir(file_1.parent) == ["file_1"]
    assert os.listdir(file_2.parent) == ["file_2"]
    assert file_1.read_text() == "2"
    assert file_2.read_text() == "2"


def test_copy_file__copy_to_different_folder(tmp_path: Path) -> None:
    (tmp_path / "file_1").write_text("3")

    with SetWorkingDirectory(tmp_path):
        copy_file("file_1", "dst/file_1")

    assert set(os.listdir(tmp_path)) == set(["file_1", "dst"])
    assert os.listdir(tmp_path / "dst") == ["file_1"]
    assert (tmp_path / "file_1").read_text() == "3"
    assert (tmp_path / "dst" / "file_1").read_text() == "3"


def test_copy_file__overwrite(tmp_path: Path) -> None:
    file_1 = tmp_path / "file_1"
    file_1.write_text("4")
    file_2 = tmp_path / "file_2"
    file_2.write_text("5")

    with SetWorkingDirectory(tmp_path):
        copy_file("file_1", "file_2")

    assert set(os.listdir(tmp_path)) == set(["file_1", "file_2"])
    assert file_1.read_text() == "4"
    assert file_2.read_text() == "4"


def test_copy_file__enoent_reraised_if_not_due_to_missing_folder() -> None:
    with pytest.raises(IOError):
        copy_file("", "./dst")


###############################################################################
###############################################################################
# Tests for 'open'


def test_open_ro__uncompressed() -> None:
    with open_ro(test_file("fasta_file.fasta")) as handle:
        assert handle.read() == ">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n"


def test_open_ro__gz() -> None:
    with open_ro(test_file("fasta_file.fasta.gz")) as handle:
        assert (
            handle.read()
            == ">This_is_GZipped_FASTA!\nACGTN\n>This_is_ALSO_GZipped_FASTA!\nCGTNA\n"
        )


def test_open_ro__bz2() -> None:
    with open_ro(test_file("fasta_file.fasta.bz2")) as handle:
        assert (
            handle.read()
            == ">This_is_BZ_FASTA!\nCGTNA\n>This_is_ALSO_BZ_FASTA!\nACGTN\n"
        )


def test_open_ro__binary__uncompressed() -> None:
    with open_ro(test_file("fasta_file.fasta"), "rb") as handle:
        assert handle.read() == b">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n"


def test_open_ro__binary__gz() -> None:
    with open_ro(test_file("fasta_file.fasta.gz"), "rb") as handle:
        assert (
            handle.read()
            == b">This_is_GZipped_FASTA!\nACGTN\n>This_is_ALSO_GZipped_FASTA!\nCGTNA\n"
        )


def test_open_ro__binary__bz2() -> None:
    with open_ro(test_file("fasta_file.fasta.bz2"), "rb") as handle:
        assert (
            handle.read()
            == b">This_is_BZ_FASTA!\nCGTNA\n>This_is_ALSO_BZ_FASTA!\nACGTN\n"
        )


class OddException(RuntimeError):
    pass


def test_open_ro__close_handle_on_error() -> None:
    mocks = Mock()
    mocks.file.read.side_effect = OddException("ARGH!")
    mocks.file.__enter__ = Mock(return_value=mocks.file)
    mocks.file.__exit__ = Mock(return_value=None)
    mocks.open.return_value = mocks.file

    with patch("builtins.open", mocks.open):
        with pytest.raises(OddException):
            open_ro("/var/abc")

    mocks.assert_has_calls(
        [
            call.open("/var/abc", "rb"),
            call.file.__enter__(),
            call.file.read(2),
            call.file.__exit__(ANY, ANY, ANY),
        ]
    )


###############################################################################
###############################################################################
# Tests for 'try_remove'


def test_try_remove(tmp_path: Path) -> None:
    fpath = tmp_path / "test.txt"
    fpath.write_text("1 2 3")
    assert try_remove(fpath)
    assert not fpath.exists()


def test_try_remove__missing(tmp_path: Path) -> None:
    fpath = tmp_path / "test.txt"
    assert not try_remove(fpath)
    assert not fpath.exists()


def test_try_remove__non_file(tmp_path: Path) -> None:
    with pytest.raises(OSError):
        try_remove(tmp_path)


###############################################################################
###############################################################################
# Tests for 'try_rmtree'


def test_try_rmtree(tmp_path: Path) -> None:
    fpath = tmp_path / "testdir"
    os.mkdir(fpath)
    (fpath / "file").write_text("1 2 3")
    assert try_rmtree(fpath)
    assert not fpath.exists()


def test_try_treedir__missing(tmp_path: Path) -> None:
    fpath = tmp_path / "testdir"
    assert not try_rmtree(fpath)
    assert not fpath.exists()


###############################################################################
###############################################################################
# Tests for 'describe_files'


def test_describe_files__no_files() -> None:
    assert describe_files(()) == "No files"


def test_describe_files__single_file() -> None:
    fpath = "/var/foo/bar"
    assert describe_files((fpath,)) == repr(fpath)


def test_describe_files__same_path_abs__3_differences() -> None:
    fpaths = ("/var/foo/bar", "/var/foo/foo")
    assert describe_files(fpaths) == "2 files in '/var/foo'"


def test_describe_files__same_path_abs__2_differences() -> None:
    fpaths = ("/var/foo/faz", "/var/foo/foo")
    assert describe_files(fpaths) == "'/var/foo/f??'"


def test_describe_files__same_path_abs__1_differences() -> None:
    fpaths = ("/var/foo/faz", "/var/foo/fao")
    assert describe_files(fpaths) == "'/var/foo/fa?'"


def test_describe_files__different_paths_abs() -> None:
    fpaths = ("/var/foo/bar", "/var/bar/foo")
    assert describe_files(fpaths) == "2 files"


def test_describe_files__same_path_rel() -> None:
    fpaths = ("var/foo/bar", "var/foo/foo")
    assert describe_files(fpaths) == "2 files in 'var/foo'"


def test_describe_files__different_paths_rel() -> None:
    fpaths = ("var/foo/bar", "var/bar/foo")
    assert describe_files(fpaths) == "2 files"


def test_describe_files__iterable() -> None:
    fpaths = iter(("/var/foo/bar", "/var/foo/foo"))
    assert describe_files(fpaths) == "2 files in '/var/foo'"


###############################################################################
###############################################################################
# Tests for 'describe_paired_files'


def test_describe_paired_files__single_file() -> None:
    fpath = "/var/foo/bar"
    assert describe_paired_files((fpath,), ()) == repr(fpath)


def test_describe_paired_files__identical_files() -> None:
    fpath = "/var/foo/bar"
    ftuple = (fpath,)
    assert describe_paired_files(ftuple, ftuple) == repr(fpath)


def test_describe_paired_files__same_path__similar_files() -> None:
    files_1 = ("foo/1_abc", "foo/1_def")
    files_2 = ("foo/1_ghi", "foo/1_jkl")
    expected = "'foo/1_???'"
    result = describe_paired_files(files_1, files_2)
    assert result == expected


def test_describe_paired_files__same_path__similar_files__different_prefixes() -> None:
    files_1 = ("foo/1_abc", "foo/1_def")
    files_2 = ("foo/2_ghi", "foo/2_jkl")
    expected = "'foo/[12]_???'"
    result = describe_paired_files(files_1, files_2)
    assert result == expected


def test_describe_paired_files__same_path__similar_files__too_different() -> None:
    files_1 = ("foo/1a_abc", "foo/1a_def")
    files_2 = ("foo/2b_ghi", "foo/2b_jkl")
    expected = "2 pair(s) of files in 'foo'"
    result = describe_paired_files(files_1, files_2)
    assert result == expected


def test_describe_paired_files__same_path__different_files() -> None:
    files_1 = ("foo/1_abc", "foo/2_def")
    files_2 = ("foo/3_ghi", "foo/4_jkl")
    expected = "2 pair(s) of files in 'foo'"
    result = describe_paired_files(files_1, files_2)
    assert result == expected


def test_describe_paired_files__same_path__different_file_lens() -> None:
    files_1 = ("foo/1_a", "foo/2_de")
    files_2 = ("foo/3_g", "foo/4_jk")
    expected = "2 pair(s) of files in 'foo'"
    result = describe_paired_files(files_1, files_2)
    assert result == expected


def test_describe_paired_files__different_path_and_files() -> None:
    files_1 = ("foo/1_abc", "bar/2_def")
    files_2 = ("zed/3_ghi", "not/4_jkl")
    expected = "2 pair(s) of files"
    result = describe_paired_files(files_1, files_2)
    assert result == expected


def test_describe_paired_files__files_1_longer() -> None:
    with pytest.raises(ValueError):
        describe_paired_files(("a", "b"), ("c",))


def test_describe_paired_files__files_2_longer() -> None:
    with pytest.raises(ValueError):
        describe_paired_files(("a",), ("b", "c"))
