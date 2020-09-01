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
from unittest.mock import call, Mock

import pytest

import paleomix.atomiccmd.pprint
from paleomix.atomiccmd.command import AtomicCmd, CmdError
from paleomix.atomiccmd.sets import ParallelCmds, SequentialCmds

_SET_CLASSES = (ParallelCmds, SequentialCmds)

###############################################################################
###############################################################################
# Properties with same expected behavior for both Parallel/SequentialCmds


@pytest.mark.parametrize("cls", _SET_CLASSES)
def test_atomicsets__properties(cls):
    cmd_mock_1 = AtomicCmd(
        ("true",),
        CHECK_A=id,
        EXEC_1="false",
        IN_1="/foo/bar/in_1.file",
        IN_2="/foo/bar/in_2.file",
        OUT_1="/bar/foo/out",
        TEMP_OUT_1="out.log",
        AUX_A="/aux/fA",
        AUX_B="/aux/fB",
    )
    cmd_mock_2 = AtomicCmd(
        ("false",),
        CHECK_A=list,
        EXEC_1="echo",
        EXEC_2="java",
        IN_1="/foo/bar/in.file",
        OUT_1="out.txt",
    )

    obj = cls([cmd_mock_1, cmd_mock_2])
    assert obj.executables == cmd_mock_1.executables | cmd_mock_2.executables
    assert obj.requirements == cmd_mock_1.requirements | cmd_mock_2.requirements
    assert obj.input_files == cmd_mock_1.input_files | cmd_mock_2.input_files
    assert obj.output_files == cmd_mock_1.output_files | cmd_mock_2.output_files
    assert (
        obj.auxiliary_files == cmd_mock_1.auxiliary_files | cmd_mock_2.auxiliary_files
    )
    assert obj.expected_temp_files == frozenset(["out", "out.txt"])
    assert (
        obj.optional_temp_files
        == cmd_mock_1.optional_temp_files | cmd_mock_2.optional_temp_files
    )


_NO_CLOBBERING_KWARGS = (
    ({"OUT_A": "/foo/out.txt"}, {"OUT_B": "/bar/out.txt"}),
    ({"OUT_A": "/foo/out.txt"}, {"TEMP_OUT_B": "out.txt"}),
    ({"OUT_A": "/foo/out.txt"}, {"OUT_STDOUT": "/bar/out.txt"}),
    ({"OUT_A": "/foo/out.txt"}, {"TEMP_OUT_STDOUT": "out.txt"}),
    ({"OUT_A": "/foo/out.txt"}, {"OUT_STDERR": "/bar/out.txt"}),
    ({"OUT_A": "/foo/out.txt"}, {"TEMP_OUT_STDERR": "out.txt"}),
)

# Ensure that commands in a set doesn't clobber eachothers OUT files
@pytest.mark.parametrize("cls", _SET_CLASSES)
@pytest.mark.parametrize("kwargs_1, kwargs_2", _NO_CLOBBERING_KWARGS)
def test_atomicsets__no_clobbering(cls, kwargs_1, kwargs_2):
    cmd_1 = AtomicCmd("true", **kwargs_1)
    cmd_2 = AtomicCmd("true", **kwargs_2)
    with pytest.raises(CmdError):
        cls([cmd_1, cmd_2])


###############################################################################
###############################################################################
# Functions with same expected behavior for both Parallel/SequentialCmds


@pytest.mark.parametrize("cls", _SET_CLASSES)
def test_atomicsets__commit(cls):
    mock = Mock()
    cmd_1 = AtomicCmd(["ls"])
    cmd_1.commit = mock.commit_1
    cmd_2 = AtomicCmd(["ls"])
    cmd_2.commit = mock.commit_2
    cmd_3 = AtomicCmd(["ls"])
    cmd_3.commit = mock.commit_3

    cls((cmd_1, cmd_2, cmd_3)).commit("xTMPx")

    assert mock.mock_calls == [
        call.commit_1("xTMPx"),
        call.commit_2("xTMPx"),
        call.commit_3("xTMPx"),
    ]


@pytest.mark.parametrize("cls", _SET_CLASSES)
def test_atomicsets__commit__remove_files_on_failure(tmp_path, cls):
    (tmp_path / "tmp").mkdir()
    out_path = tmp_path / "out"

    cmd_1 = AtomicCmd(["touch", "%(OUT_FILE)s"], OUT_FILE=str(out_path / "file1"))
    cmd_2 = AtomicCmd(["touch", "%(OUT_FILE)s"], OUT_FILE=str(out_path / "file2"))
    cmd_2.commit = Mock()
    cmd_2.commit.side_effect = OSError()

    cmdset = cls((cmd_1, cmd_2))
    cmdset.run(str(tmp_path / "tmp"))
    assert cmdset.join() == [0, 0]

    with pytest.raises(OSError):
        cmdset.commit(str(tmp_path / "tmp"))

    tmp_files = [it.name for it in (tmp_path / "tmp").iterdir()]
    assert "file1" not in tmp_files
    assert "file2" in tmp_files

    assert list((tmp_path / "out").iterdir()) == []


@pytest.mark.parametrize("cls", _SET_CLASSES)
def test_atomicsets__stdout(cls):
    cmds = cls([AtomicCmd("ls")])
    with pytest.raises(CmdError):
        cmds.stdout


@pytest.mark.parametrize("cls", _SET_CLASSES)
def test_atomicsets__terminate(cls):
    mock = Mock()
    cmd_1 = AtomicCmd(["ls"])
    cmd_1.terminate = mock.terminate_1
    cmd_2 = AtomicCmd(["ls"])
    cmd_2.terminate = mock.terminate_2
    cmd_3 = AtomicCmd(["ls"])
    cmd_3.terminate = mock.terminate_3

    cmds = cls((cmd_3, cmd_2, cmd_1))
    cmds.terminate()

    assert mock.mock_calls == [
        call.terminate_3(),
        call.terminate_2(),
        call.terminate_1(),
    ]


@pytest.mark.parametrize("cls", _SET_CLASSES)
def test_atomicsets__str__(cls):
    cmds = cls([AtomicCmd("ls")])
    assert paleomix.atomiccmd.pprint.pformat(cmds) == str(cmds)


@pytest.mark.parametrize("cls", _SET_CLASSES)
def test_atomicsets__duplicate_cmds(cls):
    cmd_1 = AtomicCmd("true")
    cmd_2 = AtomicCmd("false")
    with pytest.raises(ValueError):
        cls([cmd_1, cmd_2, cmd_1])


###############################################################################
###############################################################################
# Parallel commands


def test_parallel_commands__run():
    mock = Mock()
    cmd_1 = AtomicCmd(["ls"])
    cmd_1.run = mock.run_1
    cmd_2 = AtomicCmd(["ls"])
    cmd_2.run = mock.run_2
    cmd_3 = AtomicCmd(["ls"])
    cmd_3.run = mock.run_3

    cmds = ParallelCmds((cmd_1, cmd_2, cmd_3))
    cmds.run("xTMPx")

    assert mock.mock_calls == [
        call.run_1("xTMPx"),
        call.run_2("xTMPx"),
        call.run_3("xTMPx"),
    ]


@pytest.mark.parametrize("value", (True, False))
def test_parallel_commands__ready_single(value):
    cmd = AtomicCmd(["ls"])
    cmd.ready = Mock()
    cmd.ready.return_value = value
    cmds = ParallelCmds([cmd])
    assert cmds.ready() == value

    assert cmd.ready.mock_calls


_READY_TWO_VALUES = (
    (True, True, True),
    (False, True, False),
    (True, False, False),
    (False, False, False),
)


@pytest.mark.parametrize("first, second, result", _READY_TWO_VALUES)
def test_parallel_commands__ready_two(first, second, result):
    cmd_1 = AtomicCmd(["ls"])
    cmd_1.ready = Mock()
    cmd_1.ready.return_value = first

    cmd_2 = AtomicCmd(["ls"])
    cmd_2.ready = Mock()
    cmd_2.ready.return_value = second
    cmds = ParallelCmds([cmd_1, cmd_2])
    assert cmds.ready() == result

    assert cmd_1.ready.mock_calls
    assert bool(first) == bool(cmd_2.ready.call_count)


def test_parallel_commands__join_before_run():
    mock = Mock()
    cmd_1 = AtomicCmd(["ls"])
    cmd_1.join = mock.join_1
    cmd_2 = AtomicCmd(["ls"])
    cmd_2.join = mock.join_2
    cmd_3 = AtomicCmd(["ls"])
    cmd_3.join = mock.join_3

    cmds = ParallelCmds((cmd_3, cmd_2, cmd_1))
    assert cmds.join() == [None, None, None]

    assert mock.mock_calls == []


def test_parallel_commands__join_after_run(tmp_path):
    cmds = ParallelCmds([AtomicCmd("true") for _ in range(3)])
    cmds.run(str(tmp_path))
    assert cmds.join() == [0, 0, 0]


def _setup_mocks_for_failure(*do_mocks):
    results = []
    for do_mock in do_mocks:
        if do_mock:
            mock = AtomicCmd(("sleep", 10))
        else:
            mock = AtomicCmd("false")
        results.append(mock)
    return results


def test_parallel_commands__join_failure_1(tmp_path):
    mocks = _setup_mocks_for_failure(False, True, True)
    cmds = ParallelCmds(mocks)
    cmds.run(str(tmp_path))
    assert cmds.join() == [1, "SIGTERM", "SIGTERM"]


def test_parallel_commands__join_failure_2(tmp_path):
    mocks = _setup_mocks_for_failure(True, False, True)
    cmds = ParallelCmds(mocks)
    cmds.run(str(tmp_path))
    assert cmds.join() == ["SIGTERM", 1, "SIGTERM"]


def test_parallel_commands__join_failure_3(tmp_path):
    mocks = _setup_mocks_for_failure(True, True, False)
    cmds = ParallelCmds(mocks)
    cmds.run(str(tmp_path))
    assert cmds.join() == ["SIGTERM", "SIGTERM", 1]


def test_parallel_commands__reject_sequential():
    command = AtomicCmd(["ls"])
    seqcmd = SequentialCmds([command])
    with pytest.raises(CmdError):
        ParallelCmds([seqcmd])


def test_parallel_commands__accept_parallel():
    command = AtomicCmd(["ls"])
    parcmd = ParallelCmds([command])
    ParallelCmds([parcmd])


def test_parallel_commands__reject_noncommand():
    with pytest.raises(CmdError):
        ParallelCmds([object()])


def test_parallel_commands__reject_empty_commandset():
    with pytest.raises(CmdError):
        ParallelCmds([])


###############################################################################
###############################################################################
# Sequential commands


def test_sequential_commands__atomiccmds():
    mock = Mock()
    cmd_1 = AtomicCmd(["ls"])
    cmd_1.run = mock.run_1
    cmd_1.join = mock.join_1
    cmd_1.join.return_value = [0]
    cmd_2 = AtomicCmd(["ls"])
    cmd_2.run = mock.run_2
    cmd_2.join = mock.join_2
    cmd_2.join.return_value = [0]
    cmd_3 = AtomicCmd(["ls"])
    cmd_3.run = mock.run_3
    cmd_3.join = mock.join_3
    cmd_3.join.return_value = [0]

    cmds = SequentialCmds((cmd_1, cmd_2, cmd_3))
    assert not cmds.ready()
    cmds.run("xTMPx")
    assert cmds.ready()
    assert cmds.join() == [0, 0, 0]

    assert mock.mock_calls == [
        call.run_1("xTMPx"),
        call.join_1(),
        call.run_2("xTMPx"),
        call.join_2(),
        call.run_3("xTMPx"),
        call.join_3(),
        call.join_1(),
        call.join_2(),
        call.join_3(),
    ]


def test_sequential_commands__abort_on_error_1(tmp_path):
    cmd_1 = AtomicCmd("false")
    cmd_2 = AtomicCmd(("sleep", 10))
    cmd_3 = AtomicCmd(("sleep", 10))
    cmds = SequentialCmds([cmd_1, cmd_2, cmd_3])
    cmds.run(str(tmp_path))
    assert cmds.join() == [1, None, None]


def test_sequential_commands__abort_on_error_2(tmp_path):
    cmd_1 = AtomicCmd("true")
    cmd_2 = AtomicCmd("false")
    cmd_3 = AtomicCmd(("sleep", 10))
    cmds = SequentialCmds([cmd_1, cmd_2, cmd_3])
    cmds.run(str(tmp_path))
    assert cmds.join() == [0, 1, None]


def test_sequential_commands__abort_on_error_3(tmp_path):
    cmd_1 = AtomicCmd("true")
    cmd_2 = AtomicCmd("true")
    cmd_3 = AtomicCmd("false")
    cmds = SequentialCmds([cmd_1, cmd_2, cmd_3])
    cmds.run(str(tmp_path))
    assert cmds.join() == [0, 0, 1]


def test_sequential_commands__accept_parallel():
    command = AtomicCmd(["ls"])
    parcmd = ParallelCmds([command])
    SequentialCmds([parcmd])


def test_sequential_commands__accept_sequential():
    command = AtomicCmd(["ls"])
    seqcmd = SequentialCmds([command])
    SequentialCmds([seqcmd])


def test_sequential_commands__reject_noncommand():
    with pytest.raises(CmdError):
        SequentialCmds([object()])


def test_sequential_commands__reject_empty_commandset():
    with pytest.raises(CmdError):
        SequentialCmds([])
