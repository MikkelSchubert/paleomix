#!/usr/bin/python3
import os
import signal
import weakref

from collections import OrderedDict
from unittest.mock import call, Mock, patch
from pathlib import Path

import pytest

import paleomix.atomiccmd.command2
import paleomix.common.fileutils as fileutils

from paleomix.common.versions import RequirementObj
from paleomix.atomiccmd.command import CmdError
from paleomix.atomiccmd.command2 import (
    AtomicCmd2,
    _AtomicFile,
    AuxilleryFile,
    Executable,
    InputFile,
    OutputFile,
)

########################################################################################
# AtomicFile


_ATOMICFILE_CLASSES = (
    _AtomicFile,
    InputFile,
    OutputFile,
    Executable,
    AuxilleryFile,
)

_ATOMICFILE_INVALID_VALUES = (
    None,
    1,
    set(),
    [1, 2, 3],
    {},
    frozenset(),
    range(3),
    AtomicCmd2.DEVNULL,
    AtomicCmd2.PIPE,
)


def test_atomicfile__repr__():
    assert repr(_AtomicFile("foo/bar")) == "_AtomicFile('foo/bar')"
    assert repr(AuxilleryFile("foo/bar")) == "AuxilleryFile('foo/bar')"
    assert repr(Executable("foo/bar")) == "Executable('foo/bar')"
    assert repr(InputFile("bar")) == "InputFile('bar', False)"
    assert repr(InputFile("bar", temporary=True)) == "InputFile('bar', True)"
    assert repr(OutputFile("bar")) == "OutputFile('bar', False)"
    assert repr(OutputFile("bar", temporary=True)) == "OutputFile('bar', True)"


def test_atomicfile__valid_paths():
    assert _AtomicFile("").path == ""
    assert _AtomicFile("/foo/bar").path == "/foo/bar"
    assert _AtomicFile(Path("/foo/bar")).path == "/foo/bar"


@pytest.mark.parametrize("cls", _ATOMICFILE_CLASSES)
@pytest.mark.parametrize("value", _ATOMICFILE_INVALID_VALUES)
def test_atomiccmd2__paths__invalid_values(cls, value):
    with pytest.raises(TypeError):
        cls(value)


@pytest.mark.parametrize("cls", _ATOMICFILE_CLASSES)
def test_atomiccmd2__paths__byte_path(cls):
    with pytest.raises(ValueError):
        cls(b"/byte/path")


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


@pytest.mark.parametrize("cls", (InputFile, OutputFile))
@pytest.mark.parametrize("value", _INVALID_TEMP_PATHS)
def test_atomiccmd2__paths__invalid_temp_paths(cls, value):
    cls(value)
    with pytest.raises(ValueError):
        cls(value, temporary=True)


########################################################################################
# Constructor: Executable


def test_atomiccmd2__command_str():
    cmd = AtomicCmd2("ls")
    assert cmd.executables == frozenset(["ls"])


def test_atomiccmd2__executables_empty_str():
    with pytest.raises(ValueError, match="Empty command in AtomicCmd2 constructor"):
        AtomicCmd2("")


def test_atomiccmd2__command_tuple():
    cmd = AtomicCmd2(("cd", "."))
    assert cmd.executables == frozenset(["cd"])


def test_atomiccmd2__executables_empty_tuple():
    with pytest.raises(ValueError, match="Empty command in AtomicCmd2 constructor"):
        AtomicCmd2(())


def test_atomiccmd2__executables_empty_str_in_tuple():
    with pytest.raises(ValueError, match="Empty command in AtomicCmd2 constructor"):
        AtomicCmd2((""))


def test_atomiccmd2__executable_class():
    cmd = AtomicCmd2((Executable("cd"), "."))
    assert cmd.executables == frozenset(["cd"])


########################################################################################
# Constructor: Files in commands


# Check that specified paths/etc. are available via getters
def test_atomiccmd2__paths():
    cmd = AtomicCmd2(
        [
            "ls",
            AuxilleryFile("/path/to/index"),
            InputFile("/a/b/c"),
            OutputFile("foo/bar"),
            OutputFile("xyb", temporary=True),
        ],
    )

    assert cmd.executables == frozenset(["ls"])
    assert cmd.input_files == frozenset(["/a/b/c"])
    assert cmd.output_files == frozenset(["foo/bar"])
    assert cmd.auxiliary_files == frozenset(["/path/to/index"])


@pytest.mark.parametrize("value", _ATOMICFILE_INVALID_VALUES)
def test_atomiccmd2__values_in_path(value):
    cmd = AtomicCmd2(["echo", value])
    assert cmd.to_call("/tmp") == ["echo", str(value)]


########################################################################################
# Constructor: Extra files

# Check that specified paths/etc. are available via getters
def test_atomiccmd2__extra_files():
    cmd = AtomicCmd2(
        "ls",
        extra_files=(
            InputFile("/x/y/z"),
            InputFile("tmp_in", temporary=True),
            OutputFile("/out/foo"),
            Executable("true"),
            AuxilleryFile("wat/wat"),
        ),
    )

    assert cmd.executables == frozenset(["ls", "true"])
    assert cmd.input_files == frozenset(["/x/y/z"])
    assert cmd.output_files == frozenset(["/out/foo"])
    assert cmd.auxiliary_files == frozenset(["wat/wat"])


@pytest.mark.parametrize("value", _ATOMICFILE_INVALID_VALUES)
def test_atomiccmd2__invalid_extra_paths(value):
    with pytest.raises(ValueError):
        AtomicCmd2("ls", extra_files=[value])


########################################################################################
# Constructor: expected and optional temp files


def test_atomiccmd2__expected_and_optional_temp_files():
    cmd = AtomicCmd2(
        "ls",
        stdout="/foo/bar/data.gz",
        extra_files=(
            OutputFile("tmp_out", temporary=True),
            OutputFile("/out/foo"),
        ),
    )

    assert cmd.expected_temp_files == frozenset(["data.gz", "foo"])
    assert cmd.optional_temp_files == frozenset(
        ["pipe_ls_{}.stderr".format(id(cmd)), "tmp_out"]
    )


########################################################################################
# Constructor: Overlapping files

_OVERLAPPING_OUT_FILENAMES = (
    (OutputFile("/foo/bar/output"), OutputFile("/var/output")),
    (OutputFile("output", temporary=True), OutputFile("/var/output")),
    (OutputFile("/foo/bar/output"), OutputFile("output", temporary=True)),
)


@pytest.mark.parametrize("file1, file2", _OVERLAPPING_OUT_FILENAMES)
def test_atomiccmd2__paths__overlapping_output_1(file1, file2):
    with pytest.raises(CmdError, match="multiple output files with name 'output'"):
        AtomicCmd2(("touch", file1, file2))


@pytest.mark.parametrize("file1, file2", _OVERLAPPING_OUT_FILENAMES)
def test_atomiccmd2__paths__overlapping_output_2(file1, file2):
    with pytest.raises(CmdError, match="multiple output files with name 'output'"):
        AtomicCmd2(("touch", file1), stdout=file2)


@pytest.mark.parametrize("file1, file2", _OVERLAPPING_OUT_FILENAMES)
def test_atomiccmd2__paths__overlapping_output_3(file1, file2):
    with pytest.raises(CmdError, match="multiple output files with name 'output'"):
        AtomicCmd2(("touch"), stderr=file1, stdout=file2)


########################################################################################
# Constructor: Requirements


def test_atomicmcd__requirements():
    # RequirementObjs are the standard way to do tests
    reqobj = RequirementObj(call=("echo", "version"), search="version", checks=str)
    cmd = AtomicCmd2("true", requirements=[reqobj])
    assert cmd.requirements == frozenset([reqobj])


def test_atomicmcd__callable_as_requirements():
    cmd = AtomicCmd2("true", requirements=[bool])
    assert cmd.requirements == frozenset([bool])


def test_atomiccmd2__invalid_requirements():
    with pytest.raises(TypeError, match="requirement must be callable"):
        AtomicCmd2("ls", requirements=["ls"])


########################################################################################
# STDIN


def test_atomiccmd2__stdin_valid_values():
    AtomicCmd2("true")
    AtomicCmd2("true", stdin=None)
    AtomicCmd2("true", stdin="/path/to/file")
    AtomicCmd2("true", stdin=Path("/path/to/file"))
    AtomicCmd2("true", stdin=InputFile("/path/to/file"))
    AtomicCmd2("true", stdin=InputFile("file", temporary=True))
    AtomicCmd2("true", stdin=AtomicCmd2.DEVNULL)
    AtomicCmd2("true", stdin=AtomicCmd2.PIPE)


_INVALID_STDIN_VALUES = (
    b"/path/to/file",
    OutputFile("foo"),
    AuxilleryFile("/path/to/foo"),
    _AtomicFile("/foo/bar"),
)


@pytest.mark.parametrize("value", _INVALID_STDIN_VALUES)
def test_atomiccmd2__stdin_invalid_values(value):
    with pytest.raises(ValueError):
        AtomicCmd2("true", stdin=value)


def test_atomiccmd2__stdin_basic(tmp_path):
    fname = tmp_path / "input.fasta"
    fname.write_text(">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n")
    fname = str(fname)

    cmd = AtomicCmd2("cat", stdin=fname, stdout="result.txt")
    assert cmd.input_files == frozenset([fname])
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    result = (tmp_path / "result.txt").read_text()
    assert result == ">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n"


def test_atomiccmd2__stdin_from_temp_file(tmp_path):
    cmd = AtomicCmd2(
        "cat",
        stdin=InputFile("infile.fasta", temporary=True),
        stdout="result.txt",
    )
    assert cmd.input_files == frozenset()
    (tmp_path / "infile.fasta").write_text("a\nbc\nd")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    result = (tmp_path / "result.txt").read_text()
    assert result == "a\nbc\nd"


def test_atomiccmd2__stdin_implicit_dev_null(tmp_path):
    # STDIN should be implicitly set to /dev/null; deadlocks if not
    cmd = AtomicCmd2("cat")
    cmd.run(tmp_path)
    assert cmd.join() == [0]


@pytest.mark.parametrize("value", (None, AtomicCmd2.DEVNULL))
def test_atomiccmd2__stdin_explicit_dev_null(tmp_path, value):
    # STDIN should be implicitly set to /dev/null; deadlocks if not
    cmd = AtomicCmd2("cat", stdin=value)
    cmd.run(tmp_path)
    assert cmd.join() == [0]


########################################################################################
# STDOUT/STDERR


_VALID_STDOUT_STDERR_VALUES = [
    None,
    "/path/to/file",
    Path("/path/to/file"),
    OutputFile("/path/to/file"),
    OutputFile("file", temporary=True),
    AtomicCmd2.DEVNULL,
    AtomicCmd2.PIPE,
]


@pytest.mark.parametrize("value", _VALID_STDOUT_STDERR_VALUES)
def test_atomiccmd2__stdout_valid_values(value):
    AtomicCmd2("true", stdout=value)


@pytest.mark.parametrize("value", _VALID_STDOUT_STDERR_VALUES)
def test_atomiccmd2__stderr_valid_values(value):
    AtomicCmd2("true", stderr=value)


_INVALID_STDOUT_STDERR_VALUES = (
    b"/path/to/file",
    InputFile("foo"),
    AuxilleryFile("/path/to/foo"),
    _AtomicFile("/foo/bar"),
)


@pytest.mark.parametrize("value", _INVALID_STDOUT_STDERR_VALUES)
def test_atomiccmd2__stdout_invalid_values(value):
    with pytest.raises(ValueError):
        AtomicCmd2("true", stdout=value)


@pytest.mark.parametrize("value", _INVALID_STDOUT_STDERR_VALUES)
def test_atomiccmd2__stderr_invalid_values(value):
    with pytest.raises(ValueError):
        AtomicCmd2("true", stderr=value)


def test_atomiccmd2__stdout_implicit_filename(tmp_path):
    cmd = AtomicCmd2(("echo", "foo"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / ("pipe_echo_%i.stdout" % (id(cmd),))
    assert stdout_path.read_text() == "foo\n"
    stderr_path = tmp_path / ("pipe_echo_%i.stderr" % (id(cmd),))
    assert stderr_path.read_text() == ""


def test_atomiccmd2__stdout_explicit_filename(tmp_path):
    cmd = AtomicCmd2(("echo", "foo"), stdout="/path/to/my_output.txt")
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / "my_output.txt"
    assert stdout_path.read_text() == "foo\n"
    stderr_path = tmp_path / ("pipe_echo_%i.stderr" % (id(cmd),))
    assert stderr_path.read_text() == ""


def test_atomiccmd2__stderr_implicit_filename(tmp_path):
    cmd = AtomicCmd2(("bash", "-c", "echo foo > /dev/stderr"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / ("pipe_bash_%i.stdout" % (id(cmd),))
    assert stdout_path.read_text() == ""
    stderr_path = tmp_path / ("pipe_bash_%i.stderr" % (id(cmd),))
    assert stderr_path.read_text() == "foo\n"


def test_atomiccmd2__stderr_explicit_filename(tmp_path):
    cmd = AtomicCmd2(
        ("bash", "-c", "echo foo > /dev/stderr"),
        stderr="/path/to/my_output.txt",
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    stdout_path = tmp_path / ("pipe_bash_%i.stdout" % (id(cmd),))
    assert stdout_path.read_text() == ""
    stderr_path = tmp_path / "my_output.txt"
    assert stderr_path.read_text() == "foo\n"


def test_atomiccmd2__stdout_stderr_explicit_filename(tmp_path):
    cmd = AtomicCmd2(
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


def test_atomiccmd2__temp_dir_in_path(tmp_path):
    cmd = AtomicCmd2(("echo", "-n", "%(TEMP_DIR)s"), stdout=AtomicCmd2.PIPE)
    cmd.run(tmp_path)
    path = cmd._proc.stdout.read()
    assert tmp_path.samefile(path), (tmp_path, path)
    assert cmd.join() == [0]


def test_atomiccmd2__temp_dir_inside_path(tmp_path):
    cmd = AtomicCmd2(
        ("echo", "-n", "-Djava.io.tmpdir=%(TEMP_DIR)s"),
        stdout=AtomicCmd2.PIPE,
    )
    cmd.run(tmp_path)
    assert cmd._proc.stdout.read() == (b"-Djava.io.tmpdir=" + bytes(tmp_path))
    assert cmd.join() == [0]


########################################################################################
# Constructor: set_cwd


def test_atomiccmd2__default_cwd(tmp_path):
    cwd = os.getcwd()
    cmd = AtomicCmd2("pwd", stdout=AtomicCmd2.PIPE)
    cmd.run(tmp_path)
    assert cwd == os.getcwd()
    assert cmd._proc.stdout.read().decode() == cwd + "\n"
    assert cmd.join() == [0]


def test_atomiccmd2__set_cwd(tmp_path):
    cwd = os.getcwd()
    cmd = AtomicCmd2("pwd", stdout=AtomicCmd2.PIPE, set_cwd=True)
    cmd.run(tmp_path)
    assert cwd == os.getcwd()
    assert cmd._proc.stdout.read() == bytes(tmp_path) + b"\n"
    assert cmd.join() == [0]


# Full path when set_cwd is False, rel. path when True
_IN_OUT_PATHS_WITH_SET_CWD = [
    # cls, temporary, set_cwd, expected dir
    (InputFile, False, False, ""),
    (InputFile, False, True, "${CWD}"),
    (InputFile, True, False, "${TMP}"),
    (InputFile, True, True, ""),
    (OutputFile, False, False, "${TMP}"),
    (OutputFile, False, True, ""),
    (OutputFile, True, False, "${TMP}"),
    (OutputFile, True, True, ""),
]


@pytest.mark.parametrize("cls, temp, set_cwd, expected", _IN_OUT_PATHS_WITH_SET_CWD)
def test_atomiccmd2__set_cwd__temp_in_out(tmp_path, cls, temp, set_cwd, expected):
    cmd = AtomicCmd2(
        ("echo", "-n", cls("test_file", temporary=temp)),
        stdout=OutputFile("result.txt", temporary=True),
        set_cwd=set_cwd,
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    if expected == "${CWD}":
        expected = os.getcwd()
    elif expected == "${TMP}":
        expected = fileutils.fspath(tmp_path)

    result = (tmp_path / "result.txt").read_text()
    assert result == os.path.join(expected, "test_file")


########################################################################################
# Constructor: Piping commands


def test_atomiccmd2__piping(tmp_path):
    cmd_1 = AtomicCmd2(["echo", "-n", "#@!$^"], stdout=AtomicCmd2.PIPE)
    assert cmd_1.output_files == frozenset()
    cmd_2 = AtomicCmd2(["cat"], stdin=cmd_1, stdout="piped.txt")
    assert cmd_2.input_files == frozenset()
    cmd_1.run(tmp_path)
    cmd_2.run(tmp_path)
    assert cmd_1.join() == [0]
    assert cmd_2.join() == [0]
    assert (tmp_path / "piped.txt").read_text() == "#@!$^"


def test_atomiccmd2__piping_is_only_allowed_once(tmp_path):
    cmd_1 = AtomicCmd2(["echo", "-n", "foo\nbar"], stdout=AtomicCmd2.PIPE)
    cmd_2a = AtomicCmd2(["grep", "foo"], stdin=cmd_1)
    cmd_2b = AtomicCmd2(["grep", "bar"], stdin=cmd_1)
    cmd_1.run(tmp_path)
    cmd_2a.run(tmp_path)
    with pytest.raises(CmdError):
        cmd_2b.run(tmp_path)
    assert cmd_1.join() == [0]
    assert cmd_2a.join() == [0]
    assert cmd_2b.join() == [None]


########################################################################################
# run


def test_atomiccmd2__run__already_running(tmp_path):
    cmd = AtomicCmd2(("sleep", "10"))
    cmd.run(tmp_path)
    with pytest.raises(CmdError):
        cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]


def test_atomiccmd2__run__exception_on_missing_command(tmp_path):
    cmd = AtomicCmd2(("xyzabcefgh", "10"))
    with pytest.raises(CmdError):
        cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == [None]


def test_atomiccmd2__run__exception_on_missing_command__no_wrap(tmp_path):
    cmd = AtomicCmd2(("xyzabcefgh", "10"))
    with pytest.raises(CmdError, match="Error = FileNotFoundError"):
        cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == [None]


def test_atomiccmd2__run__invalid_temp(tmp_path):
    cmd = AtomicCmd2(("sleep", "10"))
    with pytest.raises(CmdError):
        cmd.run(tmp_path / "foo")
    cmd.terminate()
    assert cmd.join() == [None]


########################################################################################
# Ready


def test_atomiccmd2__ready_1(tmp_path):
    cmd = AtomicCmd2("ls")
    assert cmd.join() == [None]
    assert not cmd.ready()
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert cmd.ready()


def test_atomiccmd2__ready_2(tmp_path):
    cmd = AtomicCmd2(("sleep", "10"))
    cmd.run(tmp_path)
    assert not cmd.ready()
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]
    assert cmd.ready()


########################################################################################
# Join / wait


@pytest.mark.parametrize("call, after", (("true", 0), ("false", 1)))
def test_atomiccmd2__join(tmp_path, call, after):
    cmd = AtomicCmd2(call)
    assert cmd.join() == [None]
    cmd.run(tmp_path)
    assert cmd.join() == [after]


@pytest.mark.parametrize("call, after", (("true", 0), ("false", 1)))
def test_atomiccmd2__wait(tmp_path, call, after):
    cmd = AtomicCmd2(call)
    assert cmd.wait() is None
    cmd.run(tmp_path)
    assert cmd.wait() == after


########################################################################################
# Terminate


def test_atomiccmd2__terminate(tmp_path):
    cmd = AtomicCmd2(("sleep", "10"))
    cmd.run(tmp_path)

    with patch("os.killpg", wraps=os.killpg) as os_killpg:
        cmd.terminate()
        assert cmd.join() == ["SIGTERM"]

        assert os_killpg.mock_calls == [call(cmd._proc.pid, signal.SIGTERM)]


def test_atomiccmd2__terminate_exception(tmp_path):
    killpg = os.killpg
    cmd = AtomicCmd2(("sleep", "10"))
    cmd.run(tmp_path)

    def _killpg(pid, sig):
        killpg(pid, sig)
        raise OSError("Proccess not found")

    with patch("os.killpg", wraps=_killpg) as os_killpg:
        cmd.terminate()
        assert cmd.join() == ["SIGTERM"]

        assert os_killpg.mock_calls == [call(cmd._proc.pid, signal.SIGTERM)]


# Ensure that no OSException is raised, even if the command
# managed to finish before terminate was called
def test_atomiccmd2__terminate_race_condition(tmp_path):
    cmd = AtomicCmd2("true")
    cmd.run(tmp_path)
    while cmd._proc.poll() is None:
        pass
    cmd.terminate()
    assert cmd.join() == [0]


# Calling terminate on an already joined command is acceptable
def test_atomiccmd2__terminate_after_join(tmp_path):
    cmd = AtomicCmd2("true")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    cmd.terminate()
    assert cmd.join() == [0]


# Signals are translated into strings
def test_atomiccmd2__terminate_sigterm(tmp_path):
    cmd = AtomicCmd2(("sleep", "10"))
    cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]


def test_atomiccmd2__terminate_sigkill(tmp_path):
    cmd = AtomicCmd2(("sleep", "10"))
    cmd.run(tmp_path)
    cmd._proc.kill()
    assert cmd.join() == ["SIGKILL"]


########################################################################################
# commit


def _setup_for_commit(tmp_path, create_cmd=True):
    destination = tmp_path / "out"
    tmp_path = tmp_path / "tmp"
    destination.mkdir(parents=True)
    tmp_path.mkdir(parents=True)

    if not create_cmd:
        return destination, tmp_path

    cmd = AtomicCmd2(("touch", OutputFile(destination / "1234")))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    return destination, tmp_path, cmd


def test_atomiccmd2__commit_simple(tmp_path):
    destination, tmp_path, cmd = _setup_for_commit(tmp_path)
    cmd.commit(tmp_path)
    assert not (tmp_path / "1234").exists()
    assert (destination / "1234").exists()


def test_atomiccmd2__commit_temp_out(tmp_path):
    dest, temp = _setup_for_commit(tmp_path, create_cmd=False)
    cmd = AtomicCmd2(
        ("echo", "foo"),
        stdout=dest / "foo.txt",
        extra_files=[OutputFile("bar.txt", temporary=True)],
    )
    cmd.run(temp)
    assert cmd.join() == [0]
    (temp / "bar.txt").write_text("1 2 3")
    cmd.commit(temp)
    assert os.listdir(fileutils.fspath(temp)) == []
    assert os.listdir(fileutils.fspath(dest)) == ["foo.txt"]


def test_atomiccmd2__commit_temp_only(tmp_path):
    cmd = AtomicCmd2(("echo", "foo"), stdout=OutputFile("bar.txt", temporary=True))
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert (tmp_path / "bar.txt").exists()
    cmd.commit(tmp_path)
    assert os.listdir(fileutils.fspath(tmp_path)) == []


def test_atomiccmd2__commit_before_run():
    cmd = AtomicCmd2("true")
    with pytest.raises(CmdError):
        cmd.commit("/tmp")


def test_atomiccmd2__commit_while_running(tmp_path):
    cmd = AtomicCmd2(("sleep", "10"))
    cmd.run(tmp_path)
    with pytest.raises(CmdError):
        cmd.commit(tmp_path)
    cmd.terminate()
    cmd.join()


def test_atomiccmd2__commit_before_join(tmp_path):
    cmd = AtomicCmd2(("sleep", "0.1"))
    cmd.run(tmp_path)
    while cmd._proc.poll() is None:
        pass
    with pytest.raises(CmdError):
        cmd.commit(tmp_path)
    cmd.join()


# The temp path might differ, as long as the actual path is the same
def test_atomiccmd2__commit_temp_folder(tmp_path):
    destination, tmp_path, cmd = _setup_for_commit(tmp_path)
    cmd.commit(tmp_path.resolve())
    assert not (tmp_path / "1234").exists()
    assert (destination / "1234").exists()


def test_atomiccmd2__commit_wrong_temp_folder(tmp_path):
    destination, tmp_path, cmd = _setup_for_commit(tmp_path)
    with pytest.raises(CmdError):
        cmd.commit(destination)


def test_atomiccmd2__commit_missing_files(tmp_path):
    destination, tmp_path = _setup_for_commit(tmp_path, False)
    cmd = AtomicCmd2(
        ("touch", OutputFile(destination / "1234")),
        extra_files=[OutputFile(destination / "4567")],
    )
    cmd.run(tmp_path)
    cmd.join()
    before = frozenset(tmp_path.iterdir())
    with pytest.raises(CmdError):
        cmd.commit(tmp_path)
    assert before == frozenset(tmp_path.iterdir())


def test_atomiccmd2__commit_failure_cleanup(tmp_path):
    counter = []
    move_file = fileutils.move_file

    def _monkey_move_file(source, destination):
        if counter:
            raise OSError("ARRRGHHH!")
        counter.append(destination)

        return move_file(source, destination)

    destination, tmp_path = _setup_for_commit(tmp_path, False)
    command = AtomicCmd2(
        (
            "touch",
            OutputFile(destination / "file_1"),
            OutputFile(destination / "file_2"),
            OutputFile(destination / "file_3"),
        )
    )

    try:
        fileutils.move_file = _monkey_move_file
        command.run(tmp_path)
        assert command.join() == [0]
        with pytest.raises(OSError):
            command.commit(tmp_path)

        assert tuple(destination.iterdir()) == ()
    finally:
        fileutils.move_file = move_file


def test_atomiccmd2__commit_with_pipes(tmp_path):
    destination, tmp_path = _setup_for_commit(tmp_path, False)
    command_1 = AtomicCmd2(("echo", "Hello, World!"), stdout=AtomicCmd2.PIPE)
    command_2 = AtomicCmd2(("gzip",), stdin=command_1, stdout=(destination / "foo.gz"))

    command_1.run(tmp_path)
    command_2.run(tmp_path)

    assert command_1.join() == [0]
    assert command_2.join() == [0]

    command_1.commit(tmp_path)
    command_2.commit(tmp_path)

    assert list(destination.iterdir()) == [destination / "foo.gz"]
    assert list(tmp_path.iterdir()) == []


########################################################################################
# append


def test_atomiccmd2__append():
    cmd = AtomicCmd2("ls")

    assert cmd.input_files == frozenset()
    cmd.append(InputFile("/foo/bar"))
    assert cmd.input_files == frozenset(["/foo/bar"])
    assert cmd.to_call("/tmp/example") == ["ls", "/foo/bar"]


def test_atomiccmd2__append_multiple():
    cmd = AtomicCmd2("ls")

    assert cmd.output_files == frozenset()
    cmd.append("--file", OutputFile("/foo/bar"))
    assert cmd.output_files == frozenset(["/foo/bar"])
    assert cmd.to_call("/tmp/example") == ["ls", "--file", "/tmp/example/bar"]


def test_atomiccmd2__append_non_str():
    cmd = AtomicCmd2("ls")

    cmd.append("--option", 17)
    assert cmd.to_call("/tmp/example") == ["ls", "--option", "17"]


def test_atomiccmd2__append_overlapping_output__temp_and_non_temp():
    cmd = AtomicCmd2("touch")
    cmd.append(OutputFile("/foo/bar/target"))

    with pytest.raises(CmdError):
        cmd.append(OutputFile("target", temporary=True))


def test_atomiccmd2__append_overlapping_output__different_instances():
    cmd = AtomicCmd2("touch")
    cmd.append(OutputFile("/foo/bar/target"))

    with pytest.raises(CmdError):
        cmd.append(OutputFile("/foo/bar/target"))


def test_atomiccmd2__append_overlapping_output__same_instance():
    cmd = AtomicCmd2("touch")

    output_file = OutputFile("/foo/bar/target")
    cmd.append(output_file)
    cmd.append(output_file)


def test_atomiccmd2__append_non_atomic_file():
    cmd = AtomicCmd2("touch")
    assert cmd.input_files == frozenset()
    assert cmd.output_files == frozenset()
    cmd.append("target")
    assert cmd.input_files == frozenset()
    assert cmd.output_files == frozenset()


def test_atomiccmd2__append_to_running_command(tmp_path):
    with AtomicCmd2("true") as cmd:
        cmd.run(tmp_path)

        with pytest.raises(CmdError, match="cannot modify already started command"):
            cmd.append("foo")


########################################################################################
# append_options


def test_append_options__empty_lists():
    cmd = AtomicCmd2("touch")
    cmd.append_options({})

    assert cmd.to_call("${TEMP_DIR}") == ["touch"]


def test_append_options__single_options():
    options = OrderedDict()
    options["--foo"] = 1
    options["--bar"] = InputFile("/foo/bar")

    cmd = AtomicCmd2("touch")
    cmd.append_options(options)

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "--foo", "1", "--bar", "/foo/bar"]


def test_append_options__none_value():
    options = OrderedDict()
    options["--foo"] = None
    options["--bar"] = None

    cmd = AtomicCmd2("touch")
    cmd.append_options(options)

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "--foo", "--bar"]


def test_append_options__multiple_options():
    options = OrderedDict()
    options["--foo"] = [3, 2, 1]
    options["--bar"] = InputFile("/foo/bar")

    cmd = AtomicCmd2("touch")
    cmd.append_options(options)

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


def test_append_options__filted_options():
    cmd = AtomicCmd2("touch")
    cmd.append_options({"foo": 1, "--bar": InputFile("/foo/bar")})

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "--bar", "/foo/bar"]


def test_append_options__custom_filtering():
    cmd = AtomicCmd2("touch")
    cmd.append_options(
        {"foo": 1, "--bar": InputFile("/foo/bar")},
        pred=lambda s: not s.startswith("-"),
    )

    assert cmd.to_call("${TEMP_DIR}") == ["touch", "foo", "1"]


def test_append_options__invalid_key():
    cmd = AtomicCmd2("touch")
    with pytest.raises(ValueError, match="17"):
        cmd.append_options({17: "yes"})


@pytest.mark.parametrize("value", (object(),))
def test_append_options__invalid_values(value):
    cmd = AtomicCmd2("touch")
    with pytest.raises(ValueError, match=str(value)):
        cmd.append_options({"--foo": value})


@pytest.mark.parametrize("value", (object(),))
def test_append_options__invalid_values_in_list(value):
    cmd = AtomicCmd2("touch")
    with pytest.raises(ValueError, match=str(value)):
        cmd.append_options({"--foo": [1, 2, value]})


########################################################################################
# add_extra_files


def test_atomiccmd2__add_extra_files():
    cmd = AtomicCmd2("ls")

    assert cmd.input_files == frozenset()
    assert cmd.output_files == frozenset()
    cmd.add_extra_files([InputFile("/foo/bar"), OutputFile("zod", temporary=True)])
    assert cmd.input_files == frozenset(["/foo/bar"])
    assert cmd.output_files == frozenset()
    assert cmd.to_call("/tmp/example") == ["ls"]


@pytest.mark.parametrize("value", ("/path/to/file", 1, _AtomicFile("foo"), -1))
def test_atomiccmd2__add_extra_files_invalid_values(value):
    cmd = AtomicCmd2("ls")

    with pytest.raises(ValueError):
        cmd.add_extra_files([value])


def test_atomiccmd2__add_extra_files_overlapping_output_1():
    cmd = AtomicCmd2("touch")
    cmd.add_extra_files([OutputFile("/foo/bar/target")])

    with pytest.raises(CmdError, match="multiple output files with name 'target'"):
        cmd.add_extra_files([OutputFile("target", temporary=True)])


def test_atomiccmd2__add_extra_files_overlapping_output_2():
    cmd = AtomicCmd2("touch")

    with pytest.raises(CmdError, match="multiple output files with name 'target'"):
        cmd.add_extra_files(
            [OutputFile("/foo/bar/target"), OutputFile("target", temporary=True)]
        )


def test_atomiccmd2__add_extra_files_to_running_command(tmp_path):
    with AtomicCmd2("true") as cmd:
        cmd.run(tmp_path)

        with pytest.raises(CmdError, match="cannot modify already started command"):
            cmd.add_extra_files(OutputFile("foo"))


########################################################################################
# __str__

# Additional tests in atomicpp_test.py
def test_atomiccmd2__str__():
    cmd = AtomicCmd2(("echo", "test"), stdin="/foo/bar")
    assert paleomix.atomiccmd.pprint.pformat(cmd) == str(cmd)


########################################################################################
# Cleanup

# FIXME: Needs better tracking of procs
#        1. track for termination until joined
#        2. don't use weak references to avoid accidental leaks

# Test that the internal list of processes is kept clean of old objects
def test_atomiccmd2__cleanup_proc__commit(tmp_path):
    assert paleomix.atomiccmd.command._PROCS == set()
    cmd = AtomicCmd2("ls")
    cmd.run(tmp_path)

    ref = next(iter(paleomix.atomiccmd.command._PROCS))
    assert ref() == cmd._proc
    assert cmd.join() == [0]
    # Commit frees proc object
    cmd.commit(tmp_path)
    assert ref() is None

    assert ref not in paleomix.atomiccmd.command._PROCS


# Test that the internal list of processes is kept clean of old objects
def test_atomiccmd2__cleanup_proc__gc(tmp_path):
    assert paleomix.atomiccmd.command._PROCS == set()
    cmd = AtomicCmd2("ls")
    cmd.run(tmp_path)

    ref = next(iter(paleomix.atomiccmd.command._PROCS))
    assert ref() == cmd._proc
    assert cmd.join() == [0]
    # GC frees proc object
    cmd = None
    assert ref() is None

    assert ref not in paleomix.atomiccmd.command._PROCS


def test_atomiccmd2__cleanup_sigterm():
    procs = [lambda: Mock(pid=7913), lambda: Mock(pid=12345)]
    with patch("paleomix.atomiccmd.command._PROCS", procs):
        patches = Mock()
        with patch("os.killpg", patches.killpg):
            with patch("sys.exit", patches.exit):
                paleomix.atomiccmd.command._cleanup_children(signal.SIGTERM, None)

                assert patches.mock_calls == [
                    call.killpg(7913, signal.SIGTERM),
                    call.killpg(12345, signal.SIGTERM),
                    call.exit(-signal.SIGTERM),
                ]


def test_atomiccmd2__cleanup_sigterm__continues_on_exception():
    procs = [lambda: Mock(pid=7913), lambda: Mock(pid=12345)]
    with patch("paleomix.atomiccmd.command._PROCS", procs):
        patches = Mock()
        with patch("os.killpg", patches.killpg):
            with patch("sys.exit", patches.exit):
                patches.killpg.side_effect = [OSError("already killed"), None]

                paleomix.atomiccmd.command._cleanup_children(signal.SIGTERM, None)

                assert patches.mock_calls == [
                    call.killpg(7913, signal.SIGTERM),
                    call.killpg(12345, signal.SIGTERM),
                    call.exit(-signal.SIGTERM),
                ]


# Ensure that the cleanup function handles weakrefs that have been freed
def test_atomiccmd2__cleanup_sigterm__dead_weakrefs():
    dead_ref = weakref.ref(AtomicCmd2("ls"))
    with patch("paleomix.atomiccmd.command._PROCS", [dead_ref]):
        patches = Mock()
        with patch("os.killpg", patches.killpg):
            with patch("sys.exit", patches.exit):
                paleomix.atomiccmd.command._cleanup_children(signal.SIGTERM, None)

                assert patches.mock_calls == [call.exit(-signal.SIGTERM)]
