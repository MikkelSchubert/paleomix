#!/usr/bin/python3
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
import os
import signal
import weakref

from unittest.mock import call, Mock, patch

import pytest

import paleomix.atomiccmd.command
import paleomix.common.fileutils as fileutils

from paleomix.common.versions import RequirementObj
from paleomix.atomiccmd.command import AtomicCmd, CmdError


###############################################################################
###############################################################################
# Constructor: Command


def test_atomiccmd__command_str():
    cmd = AtomicCmd("ls")
    assert cmd.executables == frozenset(["ls"])


def test_atomiccmd__executables_empty_str():
    with pytest.raises(ValueError, match="Empty command in AtomicCmd constructor"):
        AtomicCmd("")


def test_atomiccmd__command_tuple():
    cmd = AtomicCmd(("cd", "."))
    assert cmd.executables == frozenset(["cd"])


def test_atomiccmd__executables_empty_tuple():
    with pytest.raises(ValueError, match="Empty command in AtomicCmd constructor"):
        AtomicCmd(())


def test_atomiccmd__executables_empty_str_in_tuple():
    with pytest.raises(ValueError, match="Empty command in AtomicCmd constructor"):
        AtomicCmd((""))


###############################################################################
###############################################################################
# Constructor: set_cwd


@pytest.mark.parametrize("set_cwd", (True, False))
def test_atomiccmd__set_cwd(tmp_path, set_cwd):
    cwd = os.getcwd()
    cmd = AtomicCmd(
        ("bash", "-c", "echo -n ${PWD}"), TEMP_OUT_STDOUT="result.txt", set_cwd=set_cwd
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert cwd == os.getcwd()

    expected = str(tmp_path) if set_cwd else cwd
    result = (tmp_path / "result.txt").read_text()
    assert os.path.samefile(expected, result), "%r != %r" % (expected, result)


# Full path when set_cwd is False, rel. path when True
@pytest.mark.parametrize("set_cwd", (True, False))
@pytest.mark.parametrize("key", ("TEMP_IN_FOO", "TEMP_OUT_FOO"))
def test_atomiccmd__set_cwd__temp_in_out(tmp_path, set_cwd, key):
    cmd = AtomicCmd(
        ("echo", "-n", "%%(%s)s" % (key,)),
        TEMP_OUT_STDOUT="result.txt",
        set_cwd=set_cwd,
        **{key: "test_file"},
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    expected = os.path.join("" if set_cwd else str(tmp_path), "test_file")
    result = (tmp_path / "result.txt").read_text()
    assert os.path.abspath(expected) == os.path.abspath(result)


###############################################################################
###############################################################################
# Constructor: Paths / pipes

# Check that specified paths/etc. are available via getters
def test_atomiccmd__paths():
    cmd = AtomicCmd(
        "ls",
        IN_AAA="/a/b/c",
        IN_AAB="/x/y/z",
        TEMP_IN_ABB="tmp_in",
        OUT_AAA="/out/foo",
        OUT_BBC="foo/bar",
        TEMP_OUT_A="xyb",
        EXEC_OTHER="true",
        AUX_WAT="wat/wat",
        CHECK_FUNC=bool,
        OUT_STDERR="/var/log/pipe.stderr",
        TEMP_OUT_STDOUT="pipe.stdout",
    )

    assert cmd.executables == frozenset(["ls", "true"])
    assert cmd.requirements == frozenset([bool])
    assert cmd.input_files == frozenset(["/a/b/c", "/x/y/z"])
    assert cmd.output_files == frozenset(
        ["/out/foo", "foo/bar", "/var/log/pipe.stderr"]
    )
    assert cmd.auxiliary_files == frozenset(["wat/wat"])
    assert cmd.expected_temp_files == frozenset(["foo", "bar", "pipe.stderr"])
    assert "xyb" in cmd.optional_temp_files
    assert "pipe.stdout" in cmd.optional_temp_files


def test_atomiccmd__paths_optional():
    cmd = AtomicCmd(["ls"], IN_OPTIONAL=None, OUT_OPTIONAL=None)
    assert cmd.input_files == frozenset()
    assert cmd.output_files == frozenset()


def test_atomiccmd__pipes_stdin(tmp_path):
    fname = tmp_path / "input.fasta"
    fname.write_text(">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n")
    fname = str(fname)

    cmd = AtomicCmd("cat", IN_STDIN=fname, OUT_STDOUT="result.txt")
    assert cmd.input_files == frozenset([fname])
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    result = (tmp_path / "result.txt").read_text()
    assert result == ">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n"


def test_atomiccmd__pipes_stdin__temp_file(tmp_path):
    cmd = AtomicCmd("cat", TEMP_IN_STDIN="infile.fasta", OUT_STDOUT="result.txt")
    assert cmd.input_files == frozenset()
    (tmp_path / "infile.fasta").write_text("a\nbc\nd")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    result = (tmp_path / "result.txt").read_text()
    assert result == "a\nbc\nd"


def test_atomiccmd__pipes_stdin__dev_null_implicit_1(tmp_path):
    # STDIN should be implicitly set to /dev/null; deadlocks if not
    cmd = AtomicCmd("cat")
    cmd.run(tmp_path)
    assert cmd.join() == [0]


def test_atomiccmd__pipes_stdin__dev_null_implicit_2(tmp_path):
    # STDIN should be implicitly set to /dev/null; deadlocks if not
    cmd = AtomicCmd("cat", IN_STDIN=None)
    cmd.run(tmp_path)
    assert cmd.join() == [0]


def test_atomiccmd__pipes_stdin__dev_null_explicit(tmp_path):
    # STDIN should be set to /dev/null; deadlocks if not
    cmd = AtomicCmd("cat", IN_STDIN=AtomicCmd.DEVNULL)
    cmd.run(tmp_path, wrap_errors=False)
    assert cmd.join() == [0]


# Test possible combinations of explicit / implicit saving of stdout/err
_TEST_OUTPUT_PIPES = (
    ("pipe_bash_{0}.stdout", "pipe_bash_{0}.stderr", {}),
    ("pipe_bash_{0}.stdout", "stderr.txt", {"OUT_STDERR": "stderr.txt"}),
    ("stdout.txt", "pipe_bash_{0}.stderr", {"OUT_STDOUT": "stdout.txt"}),
    (
        "stdout.txt",
        "stderr.txt",
        {"OUT_STDOUT": "stdout.txt", "OUT_STDERR": "stderr.txt"},
    ),
    (None, None, {"OUT_STDOUT": AtomicCmd.DEVNULL, "OUT_STDERR": AtomicCmd.DEVNULL}),
    (None, "pipe_bash_{0}.stderr", {"OUT_STDOUT": AtomicCmd.DEVNULL}),
    ("pipe_bash_{0}.stdout", None, {"OUT_STDERR": AtomicCmd.DEVNULL}),
)


@pytest.mark.parametrize("stdout, stderr, kwargs", _TEST_OUTPUT_PIPES)
def test_atomiccmd__pipes_out(tmp_path, stdout, stderr, kwargs):
    call = ("bash", "-c", "echo -n 'STDERR!' > /dev/stderr; echo -n 'STDOUT!';")

    cmd = AtomicCmd(call, **kwargs)
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    expected_files = []
    for (tmpl, text) in ((stdout, "STDOUT!"), (stderr, "STDERR!")):
        if tmpl is not None:
            fname = tmpl.format(id(cmd))
            result = (tmp_path / fname).read_text()
            assert result == text
            expected_files.append(fname)

    assert set(os.listdir(str(tmp_path))) == set(expected_files)


_MALFORMED_PATH_KEYS = (
    {"IN": "/var/foo"},  # Missing key-name #1
    {"IN_": "/var/foo"},  # Missing key-name #2
    {"TEMP_OUT": "/var/foo"},  # Missing key-name #3
    {"TEMP_OUT_": "/var/foo"},  # Missing key-name #4
    {"TEMP_OUX_FOO": "foo"},  # Invalid key-type #1
    {"INS_BAR": "foo"},  # Invalid key-type #2
)


@pytest.mark.parametrize("kwargs", _MALFORMED_PATH_KEYS)
def test_atomiccmd__paths__malformed_keys(kwargs):
    with pytest.raises(ValueError):
        AtomicCmd("true", **kwargs)


_INVALID_PATH_VALUES = (
    {"IN_FILE": 1},
    {"TEMP_IN_FILE": set()},
    {"OUT_FILE": [1, 2, 3]},
    {"TEMP_OUT_FILE": 1.0},
    {"IN_STDIN": {}},
    {"TEMP_IN_STDIN": frozenset()},
    {"OUT_STDOUT": 1.7},
    {"TEMP_OUT_STDOUT": ()},
    {"OUT_STDERR": range(3)},
    {"TEMP_OUT_STDERR": -1},
)


@pytest.mark.parametrize("kwargs", _INVALID_PATH_VALUES)
def test_atomiccmd__paths__invalid_values(kwargs):
    with pytest.raises(TypeError):
        AtomicCmd("true", **kwargs)


# Subpaths are not allowed for temp IN/OUT files, neither relative nor asbsolute
_INVALID_TEMP_PATHS = (
    # No relative paths
    {"TEMP_IN_FOO": "sub/infile"},
    {"TEMP_IN_STDIN": "sub/stdin"},
    {"TEMP_OUT_FOO": "sub/outfile"},
    {"TEMP_OUT_STDOUT": "sub/stdout"},
    {"TEMP_OUT_STDERR": "sub/stderr"},
    # No absolute paths
    {"TEMP_IN_FOO": "/tmp/sub/infile"},
    {"TEMP_IN_STDIN": "/dev/sub/stdin"},
    {"TEMP_OUT_FOO": "/etc/sub/outfile"},
    {"TEMP_OUT_STDOUT": "/var/sub/stdout"},
    {"TEMP_OUT_STDERR": "/home/sub/stderr"},
)


@pytest.mark.parametrize("kwargs", _INVALID_TEMP_PATHS)
def test_atomiccmd__paths__invalid_temp_paths(kwargs):
    with pytest.raises(ValueError):
        AtomicCmd("true", **kwargs)


_OVERLAPPING_OUT_FILENAMES = (
    {"OUT_FILE_1": "/foo/bar/outfile", "OUT_FILE_2": "/var/outfile"},
    {"TEMP_OUT_FILE_1": "outfile", "OUT_FILE_1": "/var/outfile"},
    {"OUT_FILE_1": "/foo/bar/outfile", "TEMP_OUT_FILE_1": "outfile"},
    {"TEMP_OUT_FILE_1": "outfile", "TEMP_OUT_FILE_2": "outfile"},
    {"OUT_FILE_1": "/foo/bar/outfile", "OUT_STDOUT": "/var/outfile"},
    {"TEMP_OUT_FILE_1": "outfile", "OUT_STDOUT": "/var/outfile"},
    {"OUT_FILE_1": "/foo/bar/outfile", "TEMP_OUT_STDOUT": "outfile"},
    {"TEMP_OUT_FILE_1": "outfile", "TEMP_OUT_STDOUT": "outfile"},
    {"OUT_FILE_1": "/foo/bar/outfile", "OUT_STDERR": "/var/outfile"},
    {"TEMP_OUT_FILE_1": "outfile", "OUT_STDERR": "/var/outfile"},
    {"OUT_FILE_1": "/foo/bar/outfile", "TEMP_OUT_STDERR": "outfile"},
    {"TEMP_OUT_FILE_1": "outfile", "TEMP_OUT_STDERR": "outfile"},
)

# All OUT_ files must be unique, including all TEMP_OUT_
@pytest.mark.parametrize("kwargs", _OVERLAPPING_OUT_FILENAMES)
def test_atomiccmd__paths__overlapping_output(kwargs):
    with pytest.raises(ValueError):
        AtomicCmd(("ls",), **kwargs)


# A pipe can be w/wo TEMP_, but not both
@pytest.mark.parametrize("key", ("IN_STDIN", "OUT_STDOUT", "OUT_STDERR"))
def test_atomiccmd__pipes__duplicates(key):
    kwargs = {"TEMP_" + key: "temp_file", key: "file"}
    with pytest.raises(CmdError):
        AtomicCmd(["ls"], **kwargs)


###############################################################################
###############################################################################
# CHECK_ / EXEC_

# RequirementObjs are the standard way to do tests
def test_atomicmcd__exec__reqobj():
    reqobj = RequirementObj(call=("echo", "version"), search="version", checks=str)
    cmd = AtomicCmd("true", CHECK_VERSION=reqobj)
    assert cmd.requirements == frozenset([reqobj])


# CHECK_ is expected to be a callable
def test_atomiccmd__checks__non_callable():
    with pytest.raises(TypeError, match="CHECK_FOO must be callable"):
        AtomicCmd("ls", CHECK_FOO="ls")


# EXEC_ is expected to be a string
@pytest.mark.parametrize("obj", (str, {}, 1, b"foo"))
def test_atomiccmd__exec__invalid(obj):
    with pytest.raises(TypeError, match="EXEC_FOO must be string, not "):
        AtomicCmd("true", EXEC_FOO=obj)


###############################################################################
###############################################################################
# AUX


@pytest.mark.parametrize("obj", (str, {}, 1, b"foo"))
def test_atomiccmd__aux__invalid(obj):
    with pytest.raises(TypeError, match="AUX_FOO must be string, not "):
        AtomicCmd("true", AUX_FOO=obj)


###############################################################################
###############################################################################
# Path components


def test_atomiccmd__paths_non_str(tmp_path):
    cmd = AtomicCmd(("touch", 1234), OUT_FOO="1234", set_cwd=True)
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert (tmp_path / "1234").exists()


def test_atomiccmd__paths_missing():
    with pytest.raises(CmdError, match="Value not specified for path"):
        AtomicCmd(("touch", "%(IN_FOO)s"))


def test_atomiccmd__paths_invalid():
    with pytest.raises(CmdError, match="incomplete format"):
        AtomicCmd(("touch", "%(IN_FOO)"), IN_FOO="abc")


def test_atomiccmd__paths__key(tmp_path):
    cmd = AtomicCmd(("echo", "-n", "%(TEMP_DIR)s"), OUT_STDOUT=AtomicCmd.PIPE)
    cmd.run(tmp_path)
    path = cmd._proc.stdout.read()
    assert tmp_path.samefile(path), (tmp_path, path)
    assert cmd.join() == [0]


###############################################################################
###############################################################################
# Constructor: Piping commands


def test_atomiccmd__piping(tmp_path):
    cmd_1 = AtomicCmd(["echo", "-n", "#@!$^"], OUT_STDOUT=AtomicCmd.PIPE)
    assert cmd_1.output_files == frozenset()
    cmd_2 = AtomicCmd(["cat"], IN_STDIN=cmd_1, OUT_STDOUT="piped.txt")
    assert cmd_2.input_files == frozenset()
    cmd_1.run(tmp_path)
    cmd_2.run(tmp_path)
    assert cmd_1.join() == [0]
    assert cmd_2.join() == [0]
    result = (tmp_path / "piped.txt").read_text()
    assert result == "#@!$^"


def test_atomiccmd__piping_temp(tmp_path):
    cmd_1 = AtomicCmd(["echo", "-n", "#@!$^"], TEMP_OUT_STDOUT=AtomicCmd.PIPE)
    assert cmd_1.output_files == frozenset()
    cmd_2 = AtomicCmd(["cat"], TEMP_IN_STDIN=cmd_1, OUT_STDOUT="piped.txt")
    assert cmd_2.input_files == frozenset()
    cmd_1.run(tmp_path)
    cmd_2.run(tmp_path)
    assert cmd_1.join() == [0]
    assert cmd_2.join() == [0]
    result = (tmp_path / "piped.txt").read_text()
    assert result == "#@!$^"


# Only STDOUT takes AtomicCmd.PIPE
_INVALID_PIPE_TARGETS = (
    "IN_FILE_1",
    "IN_STDIN",
    "OUT_FILE_2",
    "OUT_STDERR",
    "TEMP_IN_FILE_1",
    "TEMP_IN_STDIN",
    "TEMP_OUT_FILE_2",
    "TEMP_OUT_STDERR",
)


@pytest.mark.parametrize("key", _INVALID_PIPE_TARGETS)
def test_atomiccmd__piping__wrong_pipe(key):
    with pytest.raises(TypeError):
        AtomicCmd("ls", **{key: AtomicCmd.PIPE})


def test_atomiccmd__piping_is_only_allowed_once(tmp_path):
    cmd_1 = AtomicCmd(["echo", "-n", "foo\nbar"], OUT_STDOUT=AtomicCmd.PIPE)
    cmd_2a = AtomicCmd(["grep", "foo"], IN_STDIN=cmd_1)
    cmd_2b = AtomicCmd(["grep", "bar"], IN_STDIN=cmd_1)
    cmd_1.run(tmp_path)
    cmd_2a.run(tmp_path)
    with pytest.raises(CmdError):
        cmd_2b.run(tmp_path)
    assert cmd_1.join() == [0]
    assert cmd_2a.join() == [0]
    assert cmd_2b.join() == [None]


###############################################################################
###############################################################################
# run


def test_atomiccmd__run__already_running(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    with pytest.raises(CmdError):
        cmd.run(tmp_path)
    cmd.terminate()
    cmd.join()


def test_atomiccmd__run__exception_on_missing_command(tmp_path):
    cmd = AtomicCmd(("xyzabcefgh", "10"))
    with pytest.raises(CmdError):
        cmd.run(tmp_path)
    cmd.terminate()
    cmd.join()


def test_atomiccmd__run__exception_on_missing_command__no_wrap(tmp_path):
    cmd = AtomicCmd(("xyzabcefgh", "10"))
    with pytest.raises(OSError):
        cmd.run(tmp_path, wrap_errors=False)
    cmd.terminate()
    cmd.join()


def test_atomiccmd__run__invalid_temp(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    with pytest.raises(CmdError):
        cmd.run(tmp_path / "foo")
    cmd.terminate()
    cmd.join()


###############################################################################
###############################################################################
# Ready


def test_atomiccmd__ready(tmp_path):
    cmd = AtomicCmd("ls")
    assert cmd.join() == [None]
    assert not cmd.ready()
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert cmd.ready()


###############################################################################
###############################################################################
# Join / wait


@pytest.mark.parametrize("call, after", (("true", 0), ("false", 1)))
def test_atomiccmd__join(tmp_path, call, after):
    cmd = AtomicCmd(call)
    assert cmd.join() == [None]
    cmd.run(tmp_path)
    assert cmd.join() == [after]


@pytest.mark.parametrize("call, after", (("true", 0), ("false", 1)))
def test_atomiccmd__wait(tmp_path, call, after):
    cmd = AtomicCmd(call)
    assert cmd.wait() is None
    cmd.run(tmp_path)
    assert cmd.wait() == after


###############################################################################
###############################################################################
# Terminate


def test_atomiccmd__terminate(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)

    with patch("os.killpg", wraps=os.killpg) as os_killpg:
        cmd.terminate()
        assert cmd.join() == ["SIGTERM"]

        assert os_killpg.mock_calls == [call(cmd._proc.pid, signal.SIGTERM)]


def test_atomiccmd__terminate_exception(tmp_path):
    killpg = os.killpg
    cmd = AtomicCmd(("sleep", "10"))
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
def test_atomiccmd__terminate_race_condition(tmp_path):
    cmd = AtomicCmd("true")
    cmd.run(tmp_path)
    while cmd._proc.poll() is None:
        pass
    cmd.terminate()
    assert cmd.join() == [0]


# Calling terminate on an already joined command is acceptable
def test_atomiccmd__terminate_after_join(tmp_path):
    cmd = AtomicCmd("true")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    cmd.terminate()
    assert cmd.join() == [0]


# Signals are translated into strings
def test_atomiccmd__terminate_sigterm(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]


def test_atomiccmd__terminate_sigkill(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    cmd._proc.kill()
    assert cmd.join() == ["SIGKILL"]


###############################################################################
###############################################################################
# commit


def _setup_for_commit(tmp_path, create_cmd=True):
    destination = tmp_path / "out"
    tmp_path = tmp_path / "tmp"
    destination.mkdir(parents=True)
    tmp_path.mkdir(parents=True)

    if not create_cmd:
        return destination, tmp_path

    cmd = AtomicCmd(("touch", "%(OUT_FOO)s"), OUT_FOO=str(destination / "1234"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]

    return destination, tmp_path, cmd


def test_atomiccmd__commit_simple(tmp_path):
    destination, tmp_path, cmd = _setup_for_commit(tmp_path)
    cmd.commit(tmp_path)
    assert not (tmp_path / "1234").exists()
    assert (destination / "1234").exists()


def test_atomiccmd__commit_temp_out(tmp_path):
    dest, temp = _setup_for_commit(tmp_path, create_cmd=False)
    cmd = AtomicCmd(
        ("echo", "foo"), OUT_STDOUT=str(dest / "foo.txt"), TEMP_OUT_FOO="bar.txt",
    )
    cmd.run(temp)
    assert cmd.join() == [0]
    (temp / "bar.txt").write_text("1 2 3")
    cmd.commit(temp)
    assert os.listdir(str(temp)) == []
    assert os.listdir(str(dest)) == ["foo.txt"]


def test_atomiccmd__commit_temp_only(tmp_path):
    cmd = AtomicCmd(("echo", "foo"), TEMP_OUT_STDOUT="bar.txt")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert (tmp_path / "bar.txt").exists()
    cmd.commit(tmp_path)
    assert os.listdir(str(tmp_path)) == []


def test_atomiccmd__commit_before_run():
    cmd = AtomicCmd("true")
    with pytest.raises(CmdError):
        cmd.commit("/tmp")


def test_atomiccmd__commit_while_running(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    with pytest.raises(CmdError):
        cmd.commit(tmp_path)
    cmd.terminate()
    cmd.join()


def test_atomiccmd__commit_before_join(tmp_path):
    cmd = AtomicCmd(("sleep", "0.1"))
    cmd.run(tmp_path)
    while cmd._proc.poll() is None:
        pass
    with pytest.raises(CmdError):
        cmd.commit(tmp_path)
    cmd.join()


# The temp path might differ, as long as the actual path is the same
def test_atomiccmd__commit_temp_folder(tmp_path):
    destination, tmp_path, cmd = _setup_for_commit(tmp_path)
    cmd.commit(tmp_path.resolve())
    assert not (tmp_path / "1234").exists()
    assert (destination / "1234").exists()


def test_atomiccmd__commit_wrong_temp_folder(tmp_path):
    destination, tmp_path, cmd = _setup_for_commit(tmp_path)
    with pytest.raises(CmdError):
        cmd.commit(destination)


def test_atomiccmd__commit_missing_files(tmp_path):
    destination, tmp_path = _setup_for_commit(tmp_path, False)
    cmd = AtomicCmd(
        ("touch", "%(OUT_FOO)s"),
        OUT_FOO=str(destination / "1234"),
        OUT_BAR=str(destination / "4567"),
    )
    cmd.run(tmp_path)
    cmd.join()
    before = frozenset(tmp_path.iterdir())
    with pytest.raises(CmdError):
        cmd.commit(tmp_path)
    assert before == frozenset(tmp_path.iterdir())


def test_atomiccmd__commit_failure_cleanup(tmp_path):
    counter = []
    move_file = fileutils.move_file

    def _monkey_move_file(source, destination):
        if counter:
            raise OSError("ARRRGHHH!")
        counter.append(destination)

        return move_file(source, destination)

    destination, tmp_path = _setup_for_commit(tmp_path, False)
    command = AtomicCmd(
        ("touch", "%(OUT_FILE_1)s", "%(OUT_FILE_2)s", "%(OUT_FILE_3)s"),
        OUT_FILE_1=str(destination / "file_1"),
        OUT_FILE_2=str(destination / "file_2"),
        OUT_FILE_3=str(destination / "file_3"),
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


def test_atomiccmd__commit_with_pipes(tmp_path):
    destination, tmp_path = _setup_for_commit(tmp_path, False)
    command_1 = AtomicCmd(("echo", "Hello, World!"), OUT_STDOUT=AtomicCmd.PIPE)
    command_2 = AtomicCmd(
        ("gzip",), IN_STDIN=command_1, OUT_STDOUT=str(destination / "foo.gz")
    )

    command_1.run(tmp_path)
    command_2.run(tmp_path)

    assert command_1.join() == [0]
    assert command_2.join() == [0]

    command_1.commit(tmp_path)
    command_2.commit(tmp_path)

    assert list(destination.iterdir()) == [destination / "foo.gz"]
    assert list(tmp_path.iterdir()) == []


###############################################################################
###############################################################################
# __str__
# Additional tests in atomicpp_test.py


def test_atomiccmd__str__():
    cmd = AtomicCmd(("echo", "test"))
    assert paleomix.atomiccmd.pprint.pformat(cmd) == str(cmd)


###############################################################################
###############################################################################
# Cleanup

# FIXME: Needs better tracking of procs
#        1. track for termination until joined
#        2. don't use weak references to avoid accidental leaks

# Test that the internal list of processes is kept clean of old objects
def test_atomiccmd__cleanup_proc__commit(tmp_path):
    assert paleomix.atomiccmd.command._PROCS == set()
    cmd = AtomicCmd("ls")
    cmd.run(tmp_path)

    ref = next(iter(paleomix.atomiccmd.command._PROCS))
    assert ref() == cmd._proc
    assert cmd.join() == [0]
    # Commit frees proc object
    cmd.commit(tmp_path)
    assert ref() is None

    assert ref not in paleomix.atomiccmd.command._PROCS


# Test that the internal list of processes is kept clean of old objects
def test_atomiccmd__cleanup_proc__gc(tmp_path):
    assert paleomix.atomiccmd.command._PROCS == set()
    cmd = AtomicCmd("ls")
    cmd.run(tmp_path)

    ref = next(iter(paleomix.atomiccmd.command._PROCS))
    assert ref() == cmd._proc
    assert cmd.join() == [0]
    # GC frees proc object
    cmd = None
    assert ref() is None

    assert ref not in paleomix.atomiccmd.command._PROCS


def test_atomiccmd__cleanup_sigterm():
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


def test_atomiccmd__cleanup_sigterm__continues_on_exception():
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
def test_atomiccmd__cleanup_sigterm__dead_weakrefs():
    dead_ref = weakref.ref(AtomicCmd("ls"))
    with patch("paleomix.atomiccmd.command._PROCS", [dead_ref]):
        patches = Mock()
        with patch("os.killpg", patches.killpg):
            with patch("sys.exit", patches.exit):
                paleomix.atomiccmd.command._cleanup_children(signal.SIGTERM, None)

                assert patches.mock_calls == [call.exit(-signal.SIGTERM)]
