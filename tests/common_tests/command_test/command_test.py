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

import os
import re
import signal
import sys
from pathlib import Path
from typing import TYPE_CHECKING, NoReturn
from unittest.mock import call, patch

import pytest

import paleomix
import paleomix.common.command
from paleomix.common import fileutils
from paleomix.common.command import (
    AtomicCmd,
    AtomicFileTypes,
    AuxiliaryFile,
    CmdError,
    Executable,
    InputFile,
    IOFileTypes,
    OutputFile,
    TempInputFile,
    TempOutputFile,
)
from paleomix.common.versions import Requirement


def escape_match(obj: object) -> str:
    return re.escape(str(obj))


########################################################################################
# AtomicFile


_ATOMICFILE_CLASSES: tuple[type[AtomicFileTypes], ...] = (
    AuxiliaryFile,
    Executable,
    InputFile,
    OutputFile,
    TempInputFile,
    TempOutputFile,
)


_ATOMICFILE_INVALID_VALUES: tuple[object, ...] = (
    None,
    1,
    {1},
    [1, 2, 3],
    {1: 2},
    frozenset((1,)),
    range(3),
    AtomicCmd.DEVNULL,
    AtomicCmd.PIPE,
)


def test_atomicfile__repr__() -> None:
    assert repr(AuxiliaryFile("foo/bar")) == "AuxiliaryFile('foo/bar')"
    assert repr(Executable("foo/bar")) == "Executable('foo/bar')"
    assert repr(InputFile("bar")) == "InputFile('bar', False)"
    assert repr(TempInputFile("bar")) == "TempInputFile('bar', True)"
    assert repr(OutputFile("bar")) == "OutputFile('bar', False)"
    assert repr(TempOutputFile("bar")) == "TempOutputFile('bar', True)"


def test_atomicfile__valid_paths() -> None:
    assert AuxiliaryFile("").path == ""
    assert AuxiliaryFile("/foo/bar").path == "/foo/bar"
    assert AuxiliaryFile(Path("/foo/bar")).path == "/foo/bar"


@pytest.mark.parametrize("cls", _ATOMICFILE_CLASSES)
@pytest.mark.parametrize("value", _ATOMICFILE_INVALID_VALUES)
def test_atomiccmd__paths__invalid_values(
    cls: type[IOFileTypes], value: object
) -> None:
    with pytest.raises(TypeError, match="expected str, bytes or os.PathLike object"):
        cls(value)  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("cls", _ATOMICFILE_CLASSES)
def test_atomiccmd__paths__byte_path(cls: type[IOFileTypes]) -> None:
    with pytest.raises(TypeError, match="invalid path"):
        cls(b"/byte/path")  # pyright: ignore[reportArgumentType]


# Subpaths are not allowed for temp IN/OUT files, neither relative nor asbsolute
_INVALID_TEMP_PATHS = (
    # No relative paths
    "sub/infile",
    "sub/stdin",
    "sub/outfile",
    "sub/stdout",
    "sub/stderr",
    # No absolute paths
    "/tmp/sub/infile",
    "/dev/sub/stdin",
    "/etc/sub/outfile",
    "/var/sub/stdout",
    "/home/sub/stderr",
)


@pytest.mark.parametrize("cls", [InputFile, OutputFile])
@pytest.mark.parametrize("value", _INVALID_TEMP_PATHS)
def test_atomiccmd__paths__invalid_temp_paths(
    cls: type[InputFile | OutputFile], value: str
) -> None:
    cls(value)
    with pytest.raises(ValueError, match="directory component in temporary path"):
        cls(value, temporary=True)


########################################################################################
# Constructor: Executable


def test_atomiccmd__command_str() -> None:
    cmd = AtomicCmd("ls")
    assert cmd.executables == frozenset(["ls"])


def test_atomiccmd__executables_empty_str() -> None:
    with pytest.raises(ValueError, match="Empty command in AtomicCmd constructor"):
        AtomicCmd("")


def test_atomiccmd__command_tuple() -> None:
    cmd = AtomicCmd(("cd", "."))
    assert cmd.executables == frozenset(["cd"])


def test_atomiccmd__executables_empty_tuple() -> None:
    with pytest.raises(ValueError, match="Empty command in AtomicCmd constructor"):
        AtomicCmd(())


def test_atomiccmd__executables_empty_str_in_tuple() -> None:
    with pytest.raises(ValueError, match="Empty command in AtomicCmd constructor"):
        AtomicCmd("")


def test_atomiccmd__executable_class() -> None:
    cmd = AtomicCmd((Executable("cd"), "."))
    assert cmd.executables == frozenset(["cd"])
    assert cmd.to_call("%(TEMP_DIR)s") == ["cd", "."]


def test_atomiccmd__command_io_types() -> None:
    with pytest.raises(TypeError, match="exe must be str or Executable"):
        AtomicCmd([InputFile("ls")])


########################################################################################
# Constructor: Files in commands


# Check that specified paths/etc. are available via getters
def test_atomiccmd__paths() -> None:
    cmd = AtomicCmd(
        [
            "ls",
            AuxiliaryFile("/path/to/index"),
            InputFile("/a/b/c"),
            OutputFile("foo/bar"),
            TempOutputFile("xyb"),
        ],
    )

    assert cmd.executables == frozenset(["ls"])
    assert cmd.input_files == frozenset(["/a/b/c"])
    assert cmd.output_files == frozenset(["foo/bar"])
    assert cmd.auxiliary_files == frozenset(["/path/to/index"])


@pytest.mark.parametrize("value", _ATOMICFILE_INVALID_VALUES)
def test_atomiccmd__values_in_path(value: object) -> None:
    # FIXME: AtomicCmd should throw on most of these types
    cmd = AtomicCmd(["echo", value])  # pyright: ignore[reportArgumentType]
    assert cmd.to_call("/tmp") == ["echo", str(value)]


########################################################################################
# Constructor: Extra files


# Check that specified paths/etc. are available via getters
def test_atomiccmd__extra_files() -> None:
    cmd = AtomicCmd(
        "ls",
        extra_files=(
            InputFile("/x/y/z"),
            TempInputFile("tmp_in"),
            OutputFile("/out/foo"),
            Executable("true"),
            AuxiliaryFile("wat/wat"),
        ),
    )

    assert cmd.executables == frozenset(["ls", "true"])
    assert cmd.input_files == frozenset(["/x/y/z"])
    assert cmd.output_files == frozenset(["/out/foo"])
    assert cmd.auxiliary_files == frozenset(["wat/wat"])


@pytest.mark.parametrize("value", _ATOMICFILE_INVALID_VALUES)
def test_atomiccmd__invalid_extra_paths(value: object) -> None:
    with pytest.raises(TypeError, match=escape_match(value)):
        AtomicCmd("ls", extra_files=[value])  # pyright: ignore[reportArgumentType]


########################################################################################
# Constructor: expected and optional temp files


def test_atomiccmd__expected_and_optional_temp_files() -> None:
    cmd = AtomicCmd(
        "ls",
        stdout="/foo/bar/data.gz",
        extra_files=(
            TempOutputFile("tmp_out"),
            OutputFile("/out/foo"),
        ),
    )

    assert cmd.expected_temp_files == frozenset(["data.gz", "foo"])
    assert cmd.optional_temp_files == frozenset(
        [f"pipe_ls_{id(cmd)}.stderr", "tmp_out"]
    )


########################################################################################
# Constructor: Overlapping files

_OVERLAPPING_OUT_FILENAMES = (
    (OutputFile("/foo/bar/output"), OutputFile("/var/output")),
    (TempOutputFile("output"), OutputFile("/var/output")),
    (OutputFile("/foo/bar/output"), TempOutputFile("output")),
)


@pytest.mark.parametrize(("file1", "file2"), _OVERLAPPING_OUT_FILENAMES)
def test_atomiccmd__paths__overlapping_output_1(
    file1: OutputFile, file2: OutputFile
) -> None:
    with pytest.raises(CmdError, match="multiple output files with name 'output'"):
        AtomicCmd(("touch", file1, file2))


@pytest.mark.parametrize(("file1", "file2"), _OVERLAPPING_OUT_FILENAMES)
def test_atomiccmd__paths__overlapping_output_2(
    file1: OutputFile, file2: OutputFile
) -> None:
    with pytest.raises(CmdError, match="multiple output files with name 'output'"):
        AtomicCmd(("touch", file1), stdout=file2)


@pytest.mark.parametrize(("file1", "file2"), _OVERLAPPING_OUT_FILENAMES)
def test_atomiccmd__paths__overlapping_output_3(
    file1: OutputFile, file2: OutputFile
) -> None:
    with pytest.raises(CmdError, match="multiple output files with name 'output'"):
        AtomicCmd(("touch"), stderr=file1, stdout=file2)


########################################################################################
# Constructor: Requirements


def test_atomicmcd__requirements() -> None:
    # RequirementObjs are the standard way to do tests
    reqobj = Requirement(call=("echo", "version"), regexp="version")
    cmd = AtomicCmd("true", requirements=[reqobj])
    assert cmd.requirements == frozenset([reqobj])


def test_atomicmcd__callable_as_requirements() -> None:
    with pytest.raises(TypeError):
        AtomicCmd(
            "true",
            requirements=[bool],  # pyright: ignore[reportArgumentType]
        )


def test_atomiccmd__invalid_requirements() -> None:
    with pytest.raises(TypeError):
        AtomicCmd(
            "ls",
            requirements=["ls"],  # pyright: ignore[reportArgumentType]
        )


########################################################################################
# STDIN


def test_atomiccmd__stdin_valid_values() -> None:
    AtomicCmd("true")
    AtomicCmd("true", stdin=None)
    AtomicCmd("true", stdin="/path/to/file")
    AtomicCmd("true", stdin=Path("/path/to/file"))
    AtomicCmd("true", stdin=InputFile("/path/to/file"))
    AtomicCmd("true", stdin=TempInputFile("file"))
    AtomicCmd("true", stdin=AtomicCmd.DEVNULL)


_INVALID_STDIN_VALUES = (
    b"/path/to/file",
    OutputFile("foo"),
    AuxiliaryFile("/path/to/foo"),
    AuxiliaryFile("/foo/bar"),
    AtomicCmd.PIPE,
)


@pytest.mark.parametrize("value", _INVALID_STDIN_VALUES)
def test_atomiccmd__stdin_invalid_values(value: object) -> None:
    with pytest.raises(ValueError, match=escape_match(value)):
        AtomicCmd("true", stdin=value)  # pyright: ignore[reportArgumentType]


def test_atomiccmd__stdin_basic(tmp_path: Path) -> None:
    sub_path = tmp_path / "subfolder"
    sub_path.mkdir()

    fname = sub_path / "input.fasta"
    fname.write_text(">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n")
    fname = str(fname)

    cmd = AtomicCmd("cat", stdin=fname, stdout="result.txt")
    assert cmd.input_files == frozenset([fname])
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    result = (tmp_path / "result.txt").read_text()
    assert result == ">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n"


def test_atomiccmd__stdin_from_temp_file(tmp_path: Path) -> None:
    cmd = AtomicCmd(
        "cat",
        stdin=TempInputFile("infile.fasta"),
        stdout="result.txt",
    )
    assert cmd.input_files == frozenset()
    (tmp_path / "infile.fasta").write_text("a\nbc\nd")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    result = (tmp_path / "result.txt").read_text()
    assert result == "a\nbc\nd"


def test_atomiccmd__stdin_implicit_dev_null(tmp_path: Path) -> None:
    # STDIN should be implicitly set to /dev/null; deadlocks if not
    cmd = AtomicCmd("cat")
    cmd.run(tmp_path)
    assert cmd.join() == [0]


@pytest.mark.parametrize("value", [None, AtomicCmd.DEVNULL])
def test_atomiccmd__stdin_explicit_dev_null(tmp_path: Path, value: None | int) -> None:
    # STDIN should be implicitly set to /dev/null; deadlocks if not
    cmd = AtomicCmd("cat", stdin=value)
    cmd.run(tmp_path)
    assert cmd.join() == [0]


########################################################################################
# STDOUT/STDERR

if TYPE_CHECKING:
    OutPipeTypes = int | str | Path | OutputFile | None

_VALID_STDOUT_STDERR_VALUES: list[OutPipeTypes] = [
    None,
    "/path/to/file",
    Path("/path/to/file"),
    OutputFile("/path/to/file"),
    TempOutputFile("file"),
    AtomicCmd.DEVNULL,
]


@pytest.mark.parametrize("value", _VALID_STDOUT_STDERR_VALUES)
def test_atomiccmd__stdout_valid_values(value: OutPipeTypes) -> None:
    AtomicCmd("true", stdout=value)


def test_atomiccmd__stdout_valid_values__pipe() -> None:
    AtomicCmd("true", stdout=AtomicCmd.PIPE)


@pytest.mark.parametrize("value", _VALID_STDOUT_STDERR_VALUES)
def test_atomiccmd__stderr_valid_values(value: OutPipeTypes) -> None:
    AtomicCmd("true", stderr=value)


_INVALID_STDOUT_STDERR_VALUES = (
    b"/path/to/file",
    InputFile("foo"),
    AuxiliaryFile("/path/to/foo"),
    AuxiliaryFile("/foo/bar"),
    AtomicCmd("true"),
)


@pytest.mark.parametrize("value", _INVALID_STDOUT_STDERR_VALUES)
def test_atomiccmd__stdout_invalid_values(value: object) -> None:
    with pytest.raises(ValueError, match=escape_match(value)):
        AtomicCmd("true", stdout=value)  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("value", _INVALID_STDOUT_STDERR_VALUES)
def test_atomiccmd__stderr_invalid_values(value: object) -> None:
    with pytest.raises(ValueError, match=escape_match(value)):
        AtomicCmd("true", stderr=value)  # pyright: ignore[reportArgumentType]


def test_atomiccmd__stderr_invalid_values__pipe() -> None:
    with pytest.raises(ValueError, match=escape_match(AtomicCmd.PIPE)):
        AtomicCmd("true", stderr=AtomicCmd.PIPE)


@pytest.mark.parametrize("cls", [str, Executable])
@pytest.mark.parametrize("exe", ["echo", "/bin/echo"])
def test_atomiccmd__stdout_implicit_filename(
    exe: str, cls: type[str | Executable], tmp_path: Path
) -> None:
    cmd = AtomicCmd((cls(exe), "foo"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / ("pipe_echo_%i.stdout" % (id(cmd),))
    assert stdout_path.read_text() == "foo\n"
    stderr_path = tmp_path / ("pipe_echo_%i.stderr" % (id(cmd),))
    assert stderr_path.read_text() == ""


@pytest.mark.parametrize("cls", [str, Executable])
def test_atomiccmd__stdout_implicit_filename__python(
    cls: type[str | Executable], tmp_path: Path
) -> None:
    cmd = AtomicCmd((cls("%(PYTHON)s"), "-c", "print('foo')"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    executable = os.path.basename(sys.executable)

    stdout_path = tmp_path / ("pipe_%s_%i.stdout" % (executable, id(cmd)))
    assert stdout_path.read_text() == "foo\n"
    stderr_path = tmp_path / ("pipe_%s_%i.stderr" % (executable, id(cmd)))
    assert stderr_path.read_text() == ""


def test_atomiccmd__stdout_explicit_filename(tmp_path: Path) -> None:
    cmd = AtomicCmd(("echo", "foo"), stdout="/path/to/my_output.txt")
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / "my_output.txt"
    assert stdout_path.read_text() == "foo\n"
    stderr_path = tmp_path / ("pipe_echo_%i.stderr" % (id(cmd),))
    assert stderr_path.read_text() == ""


@pytest.mark.parametrize("exe", ["bash", "/bin/bash"])
def test_atomiccmd__stderr_implicit_filename(exe: str, tmp_path: Path) -> None:
    cmd = AtomicCmd((exe, "-c", "echo foo > /dev/stderr"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / ("pipe_bash_%i.stdout" % (id(cmd),))
    assert stdout_path.read_text() == ""
    stderr_path = tmp_path / ("pipe_bash_%i.stderr" % (id(cmd),))
    assert stderr_path.read_text() == "foo\n"


def test_atomiccmd__stderr_explicit_filename(tmp_path: Path) -> None:
    cmd = AtomicCmd(
        ("bash", "-c", "echo foo > /dev/stderr"),
        stderr="/path/to/my_output.txt",
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / ("pipe_bash_%i.stdout" % (id(cmd),))
    assert stdout_path.read_text() == ""
    stderr_path = tmp_path / "my_output.txt"
    assert stderr_path.read_text() == "foo\n"


def test_atomiccmd__stdout_stderr_explicit_filename(tmp_path: Path) -> None:
    cmd = AtomicCmd(
        ("bash", "-c", "echo foo; echo bar > /dev/stderr"),
        stdout="/path/to/my_output.txt",
        stderr="/path/to/my_errors.txt",
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / "my_output.txt"
    assert stdout_path.read_text() == "foo\n"
    stderr_path = tmp_path / "my_errors.txt"
    assert stderr_path.read_text() == "bar\n"


########################################################################################
# Path components


def test_atomiccmd__temp_dir_in_path(tmp_path: Path) -> None:
    cmd = AtomicCmd(("echo", "-n", "%(TEMP_DIR)s"), stdout=AtomicCmd.PIPE)
    cmd.run(tmp_path)
    stdout, stderr = cmd.communicate()
    assert stderr is None
    path = stdout.decode()
    assert tmp_path.samefile(path), (tmp_path, path)
    assert cmd.join() == [0]


def test_atomiccmd__temp_dir_inside_path(tmp_path: Path) -> None:
    cmd = AtomicCmd(
        ("echo", "-n", "-Djava.io.tmpdir=%(TEMP_DIR)s"),
        stdout=AtomicCmd.PIPE,
    )
    cmd.run(tmp_path)
    stdout, stderr = cmd.communicate()
    assert stdout == (b"-Djava.io.tmpdir=" + bytes(tmp_path))
    assert stderr is None
    assert cmd.join() == [0]


########################################################################################
# Constructor: set_cwd


def test_atomiccmd__default_cwd(tmp_path: Path) -> None:
    cwd = os.getcwd()
    cmd = AtomicCmd("pwd", stdout=AtomicCmd.PIPE)
    cmd.run(tmp_path)
    assert cwd == os.getcwd()
    stdout, stderr = cmd.communicate()
    assert stdout.decode() == cwd + "\n"
    assert stderr is None
    assert cmd.join() == [0]


def test_atomiccmd__set_cwd(tmp_path: Path) -> None:
    cwd = os.getcwd()
    cmd = AtomicCmd("pwd", stdout=AtomicCmd.PIPE, set_cwd=True)
    cmd.run(tmp_path)
    assert cwd == os.getcwd()
    stdout, stderr = cmd.communicate()
    assert stdout == bytes(tmp_path) + b"\n"
    assert stderr is None
    assert cmd.join() == [0]


# Full path when set_cwd is False, rel. path when True
_IN_OUT_PATHS_WITH_SET_CWD = [
    # cls, set_cwd, expected dir
    (InputFile, False, "path/to/filename"),
    (InputFile, True, os.path.join(os.getcwd(), "path/to/filename")),
    (TempInputFile, False, "${TMP}/filename"),
    (TempInputFile, True, "filename"),
    (OutputFile, False, "${TMP}/filename"),
    (OutputFile, True, "filename"),
    (TempOutputFile, False, "${TMP}/filename"),
    (TempOutputFile, True, "filename"),
    (AuxiliaryFile, False, "path/to/filename"),
    (AuxiliaryFile, True, os.path.join(os.getcwd(), "path/to/filename")),
]


@pytest.mark.parametrize(("cls", "set_cwd", "expected"), _IN_OUT_PATHS_WITH_SET_CWD)
def test_atomiccmd__set_cwd__temp_in_out(
    cls: type[IOFileTypes], *, set_cwd: bool, expected: str
) -> None:
    cmd = AtomicCmd(("touch", cls("path/to/filename")), set_cwd=set_cwd)

    assert cmd.to_call("${TMP}") == ["touch", expected]


_EXECUTABLES_WITH_SET_CWD = [
    # set_cwd, input, expected output
    (False, "exec", "exec"),
    (True, "exec", "exec"),
    (False, "/bin/ls", "/bin/ls"),
    (True, "/bin/ls", "/bin/ls"),
    (False, "./bin/ls", "./bin/ls"),
    (True, "./bin/ls", os.path.join(os.getcwd(), "bin/ls")),
]


@pytest.mark.parametrize(("set_cwd", "in_exec", "out_exec"), _EXECUTABLES_WITH_SET_CWD)
def test_atomiccmd__set_cwd__executables(
    *, set_cwd: bool, in_exec: str, out_exec: str
) -> None:
    cmd = AtomicCmd([Executable(in_exec)], set_cwd=set_cwd)

    assert cmd.to_call("${TMP}") == [out_exec]


@pytest.mark.parametrize(("set_cwd", "in_exec", "out_exec"), _EXECUTABLES_WITH_SET_CWD)
def test_atomiccmd__set_cwd__executable_str(
    *, set_cwd: bool, in_exec: str, out_exec: str
) -> None:
    cmd = AtomicCmd(in_exec, set_cwd=set_cwd)

    assert cmd.to_call("${TMP}") == [out_exec]


########################################################################################
# Constructor: Piping commands


def test_atomiccmd__piping(tmp_path: Path) -> None:
    cmd_1 = AtomicCmd(["echo", "-n", "#@!$^"], stdout=AtomicCmd.PIPE)
    assert cmd_1.output_files == frozenset()
    cmd_2 = AtomicCmd(["cat"], stdin=cmd_1, stdout="piped.txt")
    assert cmd_2.input_files == frozenset()
    cmd_1.run(tmp_path)
    cmd_2.run(tmp_path)
    assert cmd_1.join() == [0]
    assert cmd_2.join() == [0]
    assert (tmp_path / "piped.txt").read_text() == "#@!$^"


def test_atomiccmd__piping_is_only_allowed_once(tmp_path: Path) -> None:
    cmd_1 = AtomicCmd(["echo", "-n", "foo\nbar"], stdout=AtomicCmd.PIPE)
    cmd_2a = AtomicCmd(["grep", "foo"], stdin=cmd_1)
    cmd_2b = AtomicCmd(["grep", "bar"], stdin=cmd_1)
    cmd_1.run(tmp_path)
    cmd_2a.run(tmp_path)
    with pytest.raises(CmdError):
        cmd_2b.run(tmp_path)
    assert cmd_1.join() == [0]
    assert cmd_2a.join() == [0]
    assert cmd_2b.join() == [None]


def test_atomiccmd__piping_must_run_in_order(tmp_path: Path) -> None:
    cmd_1 = AtomicCmd(["echo", "-n", "foo\nbar"], stdout=AtomicCmd.PIPE)
    cmd_2 = AtomicCmd(["grep", "bar"], stdin=cmd_1)

    with pytest.raises(CmdError, match="attempted to pipe non-running command"):
        cmd_2.run(tmp_path)

    assert cmd_1.join() == [None]
    assert cmd_2.join() == [None]


def test_atomiccmd__piped_command_must_have_stdout_pipe(tmp_path: Path) -> None:
    cmd_1 = AtomicCmd(["echo", "-n", "foo\nbar"])
    cmd_2 = AtomicCmd(["grep", "bar"], stdin=cmd_1)

    cmd_1.run(tmp_path)
    with pytest.raises(CmdError, match="attempted to pipe from command without stdout"):
        cmd_2.run(tmp_path)

    assert cmd_1.join() == [0]
    assert cmd_2.join() == [None]


########################################################################################
# Constructor: Python interpreter


@pytest.mark.parametrize("set_cwd", [True, False])
def test_atomiccmd__python_call(*, set_cwd: bool) -> None:
    cmd = AtomicCmd("%(PYTHON)s", set_cwd=set_cwd)

    assert cmd.to_call("%(TMP_DIR)s") == [sys.executable]
    assert cmd.executables == frozenset(["%(PYTHON)s"])


@pytest.mark.parametrize("set_cwd", [True, False])
def test_atomiccmd__python_exec(*, set_cwd: bool) -> None:
    cmd = AtomicCmd(["ls", "%(PYTHON)s"], set_cwd=set_cwd)

    assert cmd.to_call("%(TMP_DIR)s") == ["ls", sys.executable]
    assert cmd.executables == frozenset(["ls", "%(PYTHON)s"])


########################################################################################
# run


def test_atomiccmd__run__already_running(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    with pytest.raises(CmdError):
        cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]


def test_atomiccmd__run__exception_on_missing_command(tmp_path: Path) -> None:
    cmd = AtomicCmd(("xyzabcefgh", "10"))
    with pytest.raises(CmdError):
        cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == [None]


def test_atomiccmd__run__exception_on_missing_command__no_wrap(tmp_path: Path) -> None:
    cmd = AtomicCmd(("xyzabcefgh", "10"))
    with pytest.raises(CmdError, match="Error = FileNotFoundError"):
        cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == [None]


def test_atomiccmd__run__invalid_temp(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "10"))
    with pytest.raises(CmdError):
        cmd.run(tmp_path / "foo")
    cmd.terminate()
    assert cmd.join() == [None]


########################################################################################
# Ready


def test_atomiccmd__ready_1(tmp_path: Path) -> None:
    cmd = AtomicCmd("ls")
    assert cmd.join() == [None]
    assert not cmd.ready()
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert cmd.ready()


def test_atomiccmd__ready_2(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    assert not cmd.ready()
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]
    assert cmd.ready()


########################################################################################
# Join / wait


@pytest.mark.parametrize(("call", "after"), [("true", 0), ("false", 1)])
def test_atomiccmd__join(tmp_path: Path, call: str, after: int) -> None:
    cmd = AtomicCmd(call)
    assert cmd.join() == [None]
    cmd.run(tmp_path)
    assert cmd.join() == [after]


@pytest.mark.parametrize(("call", "after"), [("true", 0), ("false", 1)])
def test_atomiccmd__wait(tmp_path: Path, call: str, after: int) -> None:
    cmd = AtomicCmd(call)
    assert cmd.wait() is None
    cmd.run(tmp_path)
    assert cmd.wait() == after


########################################################################################
# Terminate


def test_atomiccmd__terminate(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)

    pid = cmd.pid
    assert pid is not None

    with patch("os.killpg", wraps=os.killpg) as os_killpg:
        cmd.terminate()
        assert cmd.join() == ["SIGTERM"]

        assert os_killpg.mock_calls == [call(pid, signal.SIGTERM)]


def test_atomiccmd__terminate_exception(tmp_path: Path) -> None:
    killpg = os.killpg
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)

    pid = cmd.pid
    assert pid is not None

    def _killpg(pid: int, sig: int) -> NoReturn:
        killpg(pid, sig)
        raise OSError("Proccess not found")

    with patch("os.killpg", wraps=_killpg) as os_killpg:
        cmd.terminate()
        assert cmd.join() == ["SIGTERM"]

        assert os_killpg.mock_calls == [call(pid, signal.SIGTERM)]


# Ensure that no OSException is raised, even if the command
# managed to finish before terminate was called
def test_atomiccmd__terminate_race_condition(tmp_path: Path) -> None:
    cmd = AtomicCmd("true")
    cmd.run(tmp_path)
    stdout, stderr = cmd.communicate()
    assert stdout is None
    assert stderr is None
    cmd.terminate()
    assert cmd.join() == [0]


# Calling terminate on an already joined command is acceptable
def test_atomiccmd__terminate_after_join(tmp_path: Path) -> None:
    cmd = AtomicCmd("true")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    cmd.terminate()
    assert cmd.join() == [0]


# Signals are translated into strings
def test_atomiccmd__terminate_sigterm(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]


def test_atomiccmd__terminate_sigkill(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    cmd.terminate(signal=signal.SIGKILL)
    assert cmd.join() == ["SIGKILL"]


########################################################################################
# commit


def _setup_paths(tmp_path: Path) -> tuple[Path, Path]:
    destination = tmp_path / "out"
    tmp_path = tmp_path / "tmp"
    destination.mkdir(parents=True)
    tmp_path.mkdir(parents=True)

    return destination, tmp_path


def _setup_command(tmp_path: Path) -> tuple[Path, Path, AtomicCmd]:
    destination, tmp_path = _setup_paths(tmp_path)

    cmd = AtomicCmd(("touch", OutputFile(destination / "1234")))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    return destination, tmp_path, cmd


def test_atomiccmd__commit_simple(tmp_path: Path) -> None:
    destination, tmp_path, cmd = _setup_command(tmp_path)
    cmd.commit()
    assert not (tmp_path / "1234").exists()
    assert (destination / "1234").exists()


def test_atomiccmd__commit_temp_out(tmp_path: Path) -> None:
    dest, temp = _setup_paths(tmp_path)
    cmd = AtomicCmd(
        ("echo", "foo"),
        stdout=dest / "foo.txt",
        extra_files=[TempOutputFile("bar.txt")],
    )
    cmd.run(temp)
    assert cmd.join() == [0]
    (temp / "bar.txt").write_text("1 2 3")
    cmd.commit()
    assert os.listdir(fileutils.fspath(temp)) == []
    assert os.listdir(fileutils.fspath(dest)) == ["foo.txt"]


def test_atomiccmd__commit_temp_only(tmp_path: Path) -> None:
    cmd = AtomicCmd(("echo", "foo"), stdout=TempOutputFile("bar.txt"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert (tmp_path / "bar.txt").exists()
    cmd.commit()
    assert os.listdir(fileutils.fspath(tmp_path)) == []


def test_atomiccmd__commit_before_run() -> None:
    cmd = AtomicCmd("true")
    with pytest.raises(CmdError):
        cmd.commit()


def test_atomiccmd__commit_while_running(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    with pytest.raises(CmdError):
        cmd.commit()
    cmd.terminate()
    cmd.join()


def test_atomiccmd__commit_before_join(tmp_path: Path) -> None:
    cmd = AtomicCmd(("sleep", "0.1"))
    cmd.run(tmp_path)
    stdout, stderr = cmd.communicate()
    assert stdout is None
    assert stderr is None
    with pytest.raises(CmdError):
        cmd.commit()
    cmd.join()


def test_atomiccmd__commit_missing_files(tmp_path: Path) -> None:
    destination, tmp_path = _setup_paths(tmp_path)
    cmd = AtomicCmd(
        ("touch", OutputFile(destination / "1234")),
        extra_files=[OutputFile(destination / "4567")],
    )
    cmd.run(tmp_path)
    cmd.join()
    before = frozenset(tmp_path.iterdir())
    with pytest.raises(CmdError):
        cmd.commit()
    assert before == frozenset(tmp_path.iterdir())


def test_atomiccmd__commit_failure_cleanup(tmp_path: Path) -> None:
    counter: list[Path] = []
    move_file = fileutils.move_file

    def _monkey_move_file(source: Path, destination: Path) -> None:
        if counter:
            raise OSError("ARRRGHHH!")
        counter.append(destination)

        return move_file(source, destination)

    destination, tmp_path = _setup_paths(tmp_path)
    command = AtomicCmd(
        (
            "touch",
            OutputFile(destination / "file_1"),
            OutputFile(destination / "file_2"),
            OutputFile(destination / "file_3"),
        )
    )

    with patch("paleomix.common.fileutils.move_file", _monkey_move_file):
        command.run(tmp_path)
        assert command.join() == [0]
        with pytest.raises(OSError, match="ARRRGHHH!"):
            command.commit()

        assert tuple(destination.iterdir()) == ()


def test_atomiccmd__commit_with_pipes(tmp_path: Path) -> None:
    destination, tmp_path = _setup_paths(tmp_path)
    command_1 = AtomicCmd(("echo", "Hello, World!"), stdout=AtomicCmd.PIPE)
    command_2 = AtomicCmd(
        ("gzip",),
        stdin=command_1,
        stdout=destination / "foo.gz",
    )

    command_1.run(tmp_path)
    command_2.run(tmp_path)

    assert command_1.join() == [0]
    assert command_2.join() == [0]

    command_1.commit()
    command_2.commit()

    assert list(destination.iterdir()) == [destination / "foo.gz"]
    assert list(tmp_path.iterdir()) == []


########################################################################################
# append


def test_atomiccmd__append() -> None:
    cmd = AtomicCmd("ls")

    assert cmd.input_files == frozenset()
    cmd.append(InputFile("/foo/bar"))
    assert cmd.input_files == frozenset(["/foo/bar"])
    assert cmd.to_call("/tmp/example") == ["ls", "/foo/bar"]


def test_atomiccmd__append_multiple() -> None:
    cmd = AtomicCmd("ls")

    assert cmd.output_files == frozenset()
    cmd.append("--file", OutputFile("/foo/bar"))
    assert cmd.output_files == frozenset(["/foo/bar"])
    assert cmd.to_call("/tmp/example") == ["ls", "--file", "/tmp/example/bar"]


def test_atomiccmd__append_non_str() -> None:
    cmd = AtomicCmd("ls")

    cmd.append("--option", 17)
    assert cmd.to_call("/tmp/example") == ["ls", "--option", "17"]


def test_atomiccmd__append_overlapping_output__temp_and_non_temp() -> None:
    cmd = AtomicCmd("touch")
    cmd.append(OutputFile("/foo/bar/target"))

    with pytest.raises(CmdError):
        cmd.append(TempOutputFile("target"))


def test_atomiccmd__append_overlapping_output__different_instances() -> None:
    cmd = AtomicCmd("touch")
    cmd.append(OutputFile("/foo/bar/target"))

    with pytest.raises(CmdError):
        cmd.append(OutputFile("/foo/bar/target"))


def test_atomiccmd__append_overlapping_output__same_instance() -> None:
    cmd = AtomicCmd("touch")

    output_file = OutputFile("/foo/bar/target")
    cmd.append(output_file)
    cmd.append(output_file)


def test_atomiccmd__append_non_atomic_file() -> None:
    cmd = AtomicCmd("touch")
    assert cmd.input_files == frozenset()
    assert cmd.output_files == frozenset()
    cmd.append("target")
    assert cmd.input_files == frozenset()
    assert cmd.output_files == frozenset()


def test_atomiccmd__append_to_running_command(tmp_path: Path) -> None:
    with AtomicCmd("true") as cmd:
        cmd.run(tmp_path)

        with pytest.raises(CmdError, match="cannot modify already started command"):
            cmd.append("foo")


########################################################################################
# append_options


def test_append_options__empty_lists() -> None:
    cmd = AtomicCmd("touch")
    cmd.append_options({})

    assert cmd.to_call("${TEMP_DIR}") == ["touch"]


def test_append_options__single_options() -> None:
    cmd = AtomicCmd("touch")
    cmd.append_options({"--foo": 1, "--bar": InputFile("/foo/bar")})

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "--foo", "1", "--bar", "/foo/bar"]


def test_append_options__none_value() -> None:
    cmd = AtomicCmd("touch")
    cmd.append_options({"--foo": None, "--bar": None})

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "--foo", "--bar"]


def test_append_options__multiple_options() -> None:
    cmd = AtomicCmd("touch")
    cmd.append_options({"--foo": [3, 2, 1], "--bar": InputFile("/foo/bar")})

    assert cmd.to_call("${TEMP_DIR}") == [
        "touch",
        "--foo",
        "3",
        "--foo",
        "2",
        "--foo",
        "1",
        "--bar",
        "/foo/bar",
    ]


def test_append_options__filted_options() -> None:
    cmd = AtomicCmd("touch")
    cmd.append_options({"foo": 1, "--bar": InputFile("/foo/bar")})

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "--bar", "/foo/bar"]


def test_append_options__custom_filtering() -> None:
    cmd = AtomicCmd("touch")
    cmd.append_options(
        {"foo": 1, "--bar": InputFile("/foo/bar")},
        pred=lambda s: not s.startswith("-"),
    )

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "foo", "1"]


def test_append_options__invalid_type() -> None:
    cmd = AtomicCmd("touch")
    with pytest.raises(TypeError):
        cmd.append_options("-foo")  # pyright: ignore[reportArgumentType]


def test_append_options__invalid_key() -> None:
    cmd = AtomicCmd("touch")
    with pytest.raises(TypeError, match="17"):
        cmd.append_options({17: "yes"})  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("value", [object()])
def test_append_options__invalid_values(value: object) -> None:
    cmd = AtomicCmd("touch")
    with pytest.raises(TypeError, match=escape_match(value)):
        cmd.append_options({"--foo": value})  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("value", [object()])
def test_append_options__invalid_values_in_list(value: object) -> None:
    cmd = AtomicCmd("touch")
    with pytest.raises(TypeError, match=escape_match(value)):
        cmd.append_options(
            {"--foo": [1, 2, value]}  # pyright: ignore[reportArgumentType]
        )


########################################################################################
# add_extra_files


def test_atomiccmd__add_extra_files() -> None:
    cmd = AtomicCmd("ls")

    assert cmd.input_files == frozenset()
    assert cmd.output_files == frozenset()
    cmd.add_extra_files([InputFile("/foo/bar"), TempOutputFile("zod")])
    assert cmd.input_files == frozenset(["/foo/bar"])
    assert cmd.output_files == frozenset()
    assert cmd.to_call("/tmp/example") == ["ls"]


@pytest.mark.parametrize("value", ["/path/to/file", 1, (), -1])
def test_atomiccmd__add_extra_files_invalid_values(value: object) -> None:
    cmd = AtomicCmd("ls")

    with pytest.raises(TypeError, match=escape_match(value)):
        cmd.add_extra_files([value])  # pyright: ignore[reportArgumentType]


def test_atomiccmd__add_extra_files_overlapping_output_1() -> None:
    cmd = AtomicCmd("touch")
    cmd.add_extra_files([OutputFile("/foo/bar/target")])

    with pytest.raises(CmdError, match="multiple output files with name 'target'"):
        cmd.add_extra_files([TempOutputFile("target")])


def test_atomiccmd__add_extra_files_overlapping_output_2() -> None:
    cmd = AtomicCmd("touch")

    with pytest.raises(CmdError, match="multiple output files with name 'target'"):
        cmd.add_extra_files([OutputFile("/foo/bar/target"), TempOutputFile("target")])


def test_atomiccmd__add_extra_files_to_running_command(tmp_path: Path) -> None:
    with AtomicCmd("true") as cmd:
        cmd.run(tmp_path)

        with pytest.raises(CmdError, match="cannot modify already started command"):
            cmd.add_extra_files([OutputFile("foo")])


########################################################################################
# __str__


# Additional tests in atomicpp_test.py
def test_atomiccmd__str__() -> None:
    cmd = AtomicCmd(("echo", "test"), stdin="/foo/bar")
    assert paleomix.common.command.pformat(cmd) == str(cmd)
