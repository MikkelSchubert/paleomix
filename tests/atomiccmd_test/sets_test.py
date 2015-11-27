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
# Disable warnings for wierd function names
# pylint: disable=C0103
# Disable warnings caused by flexmock setups ("X is assigned to nothing")
# pylint: disable=W0106

import nose
import nose.tools
from nose.tools import assert_equal, assert_raises
from paleomix.common.testing import with_temp_folder

from flexmock import flexmock

import paleomix.atomiccmd.pprint
from paleomix.atomiccmd.command import AtomicCmd, CmdError
from paleomix.atomiccmd.sets import ParallelCmds, SequentialCmds


###############################################################################
###############################################################################
# Properties with same expected behavior for both Parallel/SequentialCmds

def test_atomicsets__properties():
    def _do_test(cls):
        cmd_mock_1 = AtomicCmd(("true",),
                               CHECK_A=id,
                               EXEC_1="false",
                               IN_1="/foo/bar/in_1.file",
                               IN_2="/foo/bar/in_2.file",
                               OUT_1="/bar/foo/out",
                               TEMP_OUT_1="out.log",
                               AUX_A="/aux/fA",
                               AUX_B="/aux/fB")
        cmd_mock_2 = AtomicCmd(("false",),
                               CHECK_A=list,
                               EXEC_1="echo",
                               EXEC_2="java",
                               IN_1="/foo/bar/in.file",
                               OUT_1="out.txt")

        obj = cls([cmd_mock_1, cmd_mock_2])
        assert_equal(obj.executables, cmd_mock_1.executables | cmd_mock_2.executables)
        assert_equal(obj.requirements, cmd_mock_1.requirements | cmd_mock_2.requirements)
        assert_equal(obj.input_files, cmd_mock_1.input_files | cmd_mock_2.input_files)
        assert_equal(obj.output_files, cmd_mock_1.output_files | cmd_mock_2.output_files)
        assert_equal(obj.auxiliary_files, cmd_mock_1.auxiliary_files | cmd_mock_2.auxiliary_files)
        assert_equal(obj.expected_temp_files, frozenset(["out", "out.txt"]))
        assert_equal(obj.optional_temp_files, cmd_mock_1.optional_temp_files | cmd_mock_2.optional_temp_files)

    for cls in (ParallelCmds, SequentialCmds):
        yield _do_test, cls


# Ensure that commands in a set doesn't clobber eachothers OUT files
def test_atomicsets__no_clobbering():
    def _do_test_atomicsets__no_clobbering(cls, kwargs_1, kwargs_2):
        cmd_1 = AtomicCmd("true", **kwargs_1)
        cmd_2 = AtomicCmd("true", **kwargs_2)
        assert_raises(CmdError, cls, [cmd_1, cmd_2])

    for cls in (ParallelCmds, SequentialCmds):
        yield _do_test_atomicsets__no_clobbering, cls, {"OUT_A": "/foo/out.txt"}, {"OUT_B": "/bar/out.txt"}
        yield _do_test_atomicsets__no_clobbering, cls, {"OUT_A": "/foo/out.txt"}, {"TEMP_OUT_B": "out.txt"}
        yield _do_test_atomicsets__no_clobbering, cls, {"OUT_A": "/foo/out.txt"}, {"OUT_STDOUT": "/bar/out.txt"}
        yield _do_test_atomicsets__no_clobbering, cls, {"OUT_A": "/foo/out.txt"}, {"TEMP_OUT_STDOUT": "out.txt"}
        yield _do_test_atomicsets__no_clobbering, cls, {"OUT_A": "/foo/out.txt"}, {"OUT_STDERR": "/bar/out.txt"}
        yield _do_test_atomicsets__no_clobbering, cls, {"OUT_A": "/foo/out.txt"}, {"TEMP_OUT_STDERR": "out.txt"}


###############################################################################
###############################################################################
# Functions with same expected behavior for both Parallel/SequentialCmds

def test_atomicsets__commit():
    def _do_test_atomicsets__commit(cls):
        mocks = []
        for _ in range(3):
            cmd_mock = flexmock(AtomicCmd(["ls"]))
            cmd_mock.should_receive('commit').with_args("xTMPx").once.ordered
            mocks.append(cmd_mock)

        cls(mocks).commit("xTMPx")

    yield _do_test_atomicsets__commit, ParallelCmds
    yield _do_test_atomicsets__commit, SequentialCmds


def test_atomicsets__stdout():
    @nose.tools.raises(CmdError)
    def _do_test_atomicsets__stdout(cls):
        cmds = cls([AtomicCmd("ls")])
        cmds.stdout

    yield _do_test_atomicsets__stdout, ParallelCmds
    yield _do_test_atomicsets__stdout, SequentialCmds


def test_atomicsets__terminate():
    def _do_test_atomicsets__terminate(cls):
        mocks = []
        for _ in reversed(range(3)):
            cmd_mock = flexmock(AtomicCmd("true"))
            cmd_mock.should_receive('terminate').with_args().once
            mocks.append(cmd_mock)
        cmds = cls(mocks)
        cmds.terminate()

    yield _do_test_atomicsets__terminate, ParallelCmds
    yield _do_test_atomicsets__terminate, SequentialCmds


def test_atomicsets__str__():
    def _do_test_atomicsets__str__(cls):
        cmds = cls([AtomicCmd("ls")])
        assert_equal(paleomix.atomiccmd.pprint.pformat(cmds), str(cmds))

    yield _do_test_atomicsets__str__, ParallelCmds
    yield _do_test_atomicsets__str__, SequentialCmds


def test_atomicsets__duplicate_cmds():
    def _do_test_atomicsets__duplicate_cmds(cls):
        cmd_1 = AtomicCmd("true")
        cmd_2 = AtomicCmd("false")
        assert_raises(ValueError, cls, [cmd_1, cmd_2, cmd_1])

    yield _do_test_atomicsets__duplicate_cmds, ParallelCmds
    yield _do_test_atomicsets__duplicate_cmds, SequentialCmds


###############################################################################
###############################################################################
# Parallel commands

def test_parallel_commands__run():
    mocks = []
    for _ in range(3):
        cmd_mock = flexmock(AtomicCmd(["ls"]))
        cmd_mock.should_receive('run').with_args("xTMPx").once
        mocks.append(cmd_mock)

    cmds = ParallelCmds(mocks)
    cmds.run("xTMPx")


def test_parallel_commands__ready_single():
    def _do_test_parallel_commands__ready_single(value):
        cmd_mock = flexmock(AtomicCmd(["ls"]))
        cmd_mock.should_receive('ready').and_return(value).at_least.once
        cmds = ParallelCmds([cmd_mock])
        assert_equal(cmds.ready(), value)

    yield _do_test_parallel_commands__ready_single, True
    yield _do_test_parallel_commands__ready_single, False


def test_parallel_commands__ready_two():
    def _do_test_parallel_commands__ready_two(first, second, result):
        cmd_mock_1 = flexmock(AtomicCmd(["ls"]))
        cmd_mock_1.should_receive('ready').and_return(first).at_least.once
        cmd_mock_2 = flexmock(AtomicCmd(["ls"]))
        cmd_mock_2.should_receive('ready').and_return(second)
        cmds = ParallelCmds([cmd_mock_1, cmd_mock_2])
        assert_equal(cmds.ready(), result)

    yield _do_test_parallel_commands__ready_two, True,   True,  True
    yield _do_test_parallel_commands__ready_two, False,  True, False
    yield _do_test_parallel_commands__ready_two, True,  False, False
    yield _do_test_parallel_commands__ready_two, False, False, False


def test_parallel_commands__join_before_run():
    mocks = []
    for value in reversed(range(3)):
        cmd_mock = flexmock(AtomicCmd("true"))
        cmd_mock.should_receive('join').and_return([value]).never
        mocks.append(cmd_mock)
    cmds = ParallelCmds(mocks)
    assert_equal(cmds.join(), [None, None, None])


@with_temp_folder
def test_parallel_commands__join_after_run(temp_folder):
    cmds = ParallelCmds([AtomicCmd("true") for _ in range(3)])
    cmds.run(temp_folder)
    assert_equal(cmds.join(), [0, 0, 0])


def _setup_mocks_for_failure(*do_mocks):
    results = []
    for do_mock in do_mocks:
        if do_mock:
            mock = flexmock(AtomicCmd(("sleep", 10)))
            mock.should_receive('terminate')
            mock.should_receive('join').and_return(['SIGTERM'])
        else:
            mock = AtomicCmd("false")
        results.append(mock)
    return results


@with_temp_folder
def test_parallel_commands__join_failure_1(temp_folder):
    mocks = _setup_mocks_for_failure(False, True, True)
    cmds = ParallelCmds(mocks)
    cmds.run(temp_folder)
    assert_equal(cmds.join(), [1, 'SIGTERM', 'SIGTERM'])


@with_temp_folder
def test_parallel_commands__join_failure_2(temp_folder):
    mocks = _setup_mocks_for_failure(True, False, True)
    cmds = ParallelCmds(mocks)
    cmds.run(temp_folder)
    assert_equal(cmds.join(), ['SIGTERM', 1, 'SIGTERM'])


@with_temp_folder
def test_parallel_commands__join_failure_3(temp_folder):
    mocks = _setup_mocks_for_failure(True, True, False)
    cmds = ParallelCmds(mocks)
    cmds.run(temp_folder)
    assert_equal(cmds.join(), ['SIGTERM', 'SIGTERM', 1])


def test_parallel_commands__reject_sequential():
    command = AtomicCmd(["ls"])
    seqcmd = SequentialCmds([command])
    assert_raises(CmdError, ParallelCmds, [seqcmd])


def test_parallel_commands__accept_parallel():
    command = AtomicCmd(["ls"])
    parcmd = ParallelCmds([command])
    ParallelCmds([parcmd])


@nose.tools.raises(CmdError)
def test_parallel_commands__reject_noncommand():
    ParallelCmds([object()])


@nose.tools.raises(CmdError)
def test_parallel_commands__reject_empty_commandset():
    ParallelCmds([])


###############################################################################
###############################################################################
# Sequential commands

def test_sequential_commands__atomiccmds():
    mocks = []
    for _ in range(3):
        cmd_mock = flexmock(AtomicCmd(["ls"]))
        cmd_mock.should_receive('run').with_args("xTMPx").once
        cmd_mock.should_receive('join').with_args().and_return([0]).twice
        mocks.append(cmd_mock)

    cmds = SequentialCmds(mocks)
    assert not cmds.ready()
    cmds.run("xTMPx")
    assert cmds.ready()
    assert_equal(cmds.join(), [0, 0, 0])


@with_temp_folder
@nose.tools.timed(1)
def test_sequential_commands__abort_on_error_1(temp_folder):
    cmd_1 = AtomicCmd("false")
    cmd_2 = AtomicCmd(("sleep", 10))
    cmd_3 = AtomicCmd(("sleep", 10))
    cmds = SequentialCmds([cmd_1, cmd_2, cmd_3])
    cmds.run(temp_folder)
    assert_equal(cmds.join(), [1, None, None])


@with_temp_folder
@nose.tools.timed(1)
def test_sequential_commands__abort_on_error_2(temp_folder):
    cmd_1 = AtomicCmd("true")
    cmd_2 = AtomicCmd("false")
    cmd_3 = AtomicCmd(("sleep", 10))
    cmds = SequentialCmds([cmd_1, cmd_2, cmd_3])
    cmds.run(temp_folder)
    assert_equal(cmds.join(), [0, 1, None])


@with_temp_folder
@nose.tools.timed(1)
def test_sequential_commands__abort_on_error_3(temp_folder):
    cmd_1 = AtomicCmd("true")
    cmd_2 = AtomicCmd("true")
    cmd_3 = AtomicCmd("false")
    cmds = SequentialCmds([cmd_1, cmd_2, cmd_3])
    cmds.run(temp_folder)
    assert_equal(cmds.join(), [0, 0, 1])


def test_sequential_commands__accept_parallel():
    command = AtomicCmd(["ls"])
    parcmd = ParallelCmds([command])
    SequentialCmds([parcmd])


def test_sequential_commands__accept_sequential():
    command = AtomicCmd(["ls"])
    seqcmd = SequentialCmds([command])
    SequentialCmds([seqcmd])


@nose.tools.raises(CmdError)
def test_sequential_commands__reject_noncommand():
    SequentialCmds([object()])


@nose.tools.raises(CmdError)
def test_sequential_commands__reject_empty_commandset():
    SequentialCmds([])
