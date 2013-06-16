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
# Disable warnings for wierd function names and touching private vars/funcs:
# pylint: disable=C0103
# pylint: disable=W0212

import os
import types
import signal
import weakref

import nose
from nose.tools import assert_equal, assert_in # pylint: disable=E0611
from tests.common.utils import with_temp_folder, monkeypatch, \
     get_file_contents, \
     set_file_contents

import pypeline.atomiccmd
import pypeline.atomicpp
from pypeline.common.versions import RequirementObj
from pypeline.atomiccmd import AtomicCmd, CmdError



################################################################################
################################################################################
## Constructor: Command

def test_atomiccmd__command_str():
    cmd = AtomicCmd("ls")
    assert_equal(cmd.executables, frozenset(["ls"]))

@nose.tools.raises(ValueError)
def test_atomiccmd__executables_empty_str():
    AtomicCmd("")

def test_atomiccmd__command_tuple():
    cmd = AtomicCmd(("cd", "."))
    assert_equal(cmd.executables, frozenset(["cd"]))

@nose.tools.raises(ValueError)
def test_atomiccmd__executables_empty_tuple():
    AtomicCmd(())

@nose.tools.raises(ValueError)
def test_atomiccmd__executables_empty_str_in_tuple():
    AtomicCmd((""))




################################################################################
################################################################################
## Constructor: set_cwd

def test_atomiccmd__set_cwd():
    @with_temp_folder
    def _do_test_atomiccmd__set_cwd(temp_folder, set_cwd):
        cwd = os.getcwd()
        cmd = AtomicCmd(("bash", "-c", "echo -n ${PWD}"),
                        TEMP_OUT_STDOUT = "result.txt",
                        set_cwd         = set_cwd)
        cmd.run(temp_folder)
        assert_equal(cmd.join(), [0])
        assert_equal(cwd, os.getcwd())

        expected = temp_folder if set_cwd else cwd
        result = get_file_contents(os.path.join(temp_folder, "result.txt"))
        assert os.path.samefile(expected, result), "%r != %r" % (expected, result)

    yield _do_test_atomiccmd__set_cwd, False
    yield _do_test_atomiccmd__set_cwd, True


# Full path when set_cwd is False, rel. path when True
def test_atomiccmd__set_cwd__temp_in_out():
    @with_temp_folder
    def _do_test_atomiccmd__paths_temp_in(temp_folder, set_cwd, kwargs):
        cmd = AtomicCmd(("echo", "-n", "%%(%s)s" % tuple(kwargs.keys())),
                        TEMP_OUT_STDOUT = "result.txt",
                        set_cwd         = set_cwd,
                        **kwargs)
        cmd.run(temp_folder)
        assert_equal(cmd.join(), [0])

        expected = os.path.join("" if set_cwd else temp_folder, "test_file")
        result = get_file_contents(os.path.join(temp_folder, "result.txt"))
        assert_equal(os.path.abspath(expected), os.path.abspath(result))

    yield _do_test_atomiccmd__paths_temp_in, True,  {"TEMP_IN_FOO"  : "test_file"}
    yield _do_test_atomiccmd__paths_temp_in, False, {"TEMP_IN_FOO"  : "test_file"}
    yield _do_test_atomiccmd__paths_temp_in, True,  {"TEMP_OUT_FOO" : "test_file"}
    yield _do_test_atomiccmd__paths_temp_in, False, {"TEMP_OUT_FOO" : "test_file"}



################################################################################
################################################################################
## Constructor: Paths / pipes

# Check that specified paths/etc. are available via getters
def test_atomiccmd__paths():
    cmd = AtomicCmd("ls",
                    IN_AAA      = "/a/b/c",
                    IN_AAB      = "/x/y/z",
                    TEMP_IN_ABB = "tmp_in",
                    OUT_AAA     = "/out/foo",
                    OUT_BBC     = "foo/bar",
                    TEMP_OUT_A  = "xyb",
                    EXEC_OTHER  = "true",
                    AUX_WAT     = "wat/wat",
                    CHECK_FUNC  = bool,
                    OUT_STDERR  = "/var/log/pipe.stderr",
                    TEMP_OUT_STDOUT = "pipe.stdout")

    assert_equal(cmd.executables,     frozenset(["ls", "true"]))
    assert_equal(cmd.requirements,    frozenset([bool]))
    assert_equal(cmd.input_files,     frozenset(["/a/b/c", "/x/y/z"]))
    assert_equal(cmd.output_files,    frozenset(["/out/foo", "foo/bar", "/var/log/pipe.stderr"]))
    assert_equal(cmd.auxiliary_files, frozenset(["wat/wat"]))
    assert_equal(cmd.expected_temp_files, frozenset(["foo", "bar", "pipe.stderr"]))
    assert_in("xyb", cmd.optional_temp_files)
    assert_in("pipe.stdout", cmd.optional_temp_files)


def test_atomiccmd__paths_optional():
    cmd = AtomicCmd(["ls"],
                    IN_OPTIONAL = None,
                    OUT_OPTIONAL = None)
    assert_equal(cmd.input_files, frozenset())
    assert_equal(cmd.output_files, frozenset())


@with_temp_folder
def test_atomiccmd__pipes_stdin(temp_folder):
    fname = os.path.join("tests", "data", "fasta_file.fasta")
    cmd = AtomicCmd("cat",
                    IN_STDIN   = fname,
                    OUT_STDOUT = "result.txt")
    assert_equal(cmd.input_files, frozenset([fname]))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    result = get_file_contents(os.path.join(temp_folder, "result.txt"))
    assert_equal(result, ">This_is_FASTA!\nACGTN\n>This_is_ALSO_FASTA!\nCGTNA\n")


@with_temp_folder
def test_atomiccmd__pipes_stdin__temp_file(temp_folder):
    cmd = AtomicCmd("cat",
                    TEMP_IN_STDIN  = "infile.fasta",
                    OUT_STDOUT     = "result.txt")
    assert_equal(cmd.input_files, frozenset())
    set_file_contents(os.path.join(temp_folder, "infile.fasta"), "a\nbc\nd")
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    result = get_file_contents(os.path.join(temp_folder, "infile.fasta"))
    assert_equal(result, "a\nbc\nd")


# Test possible combinations of explicit / implicit saving of stdout/err
def test_atomiccmd__pipes_out():
    @with_temp_folder
    def _do_test_atomiccmd__pipes_out(temp_folder, stdout, stderr, kwargs):
        cmd = AtomicCmd(("bash", "-c", "echo -n 'STDERR!' > /dev/stderr; echo -n 'STDOUT!';"), **kwargs)
        cmd.run(temp_folder)
        assert_equal(cmd.join(), [0])
        result_out = get_file_contents(os.path.join(temp_folder, stdout.format(id(cmd))))
        result_err = get_file_contents(os.path.join(temp_folder, stderr.format(id(cmd))))
        assert_equal(result_out, "STDOUT!")
        assert_equal(result_err, "STDERR!")

    yield _do_test_atomiccmd__pipes_out, "pipe_bash_{}.stdout", "pipe_bash_{}.stderr", {}
    yield _do_test_atomiccmd__pipes_out, "pipe_bash_{}.stdout", "stderr.txt", {"OUT_STDERR" : "stderr.txt"}
    yield _do_test_atomiccmd__pipes_out, "stdout.txt", "pipe_bash_{}.stderr", {"OUT_STDOUT" : "stdout.txt"}
    yield _do_test_atomiccmd__pipes_out, "stdout.txt", "stderr.txt", {"OUT_STDOUT" : "stdout.txt",
                                                                      "OUT_STDERR" : "stderr.txt"}

def test_atomiccmd__paths__malformed_keys():
    @nose.tools.raises(ValueError)
    def _do_test_atomiccmd__paths__malformed(kwargs):
        AtomicCmd("true", **kwargs)

    yield _do_test_atomiccmd__paths__malformed, {"IN"  : "/var/foo"} # Missing key-name #1
    yield _do_test_atomiccmd__paths__malformed, {"IN_" : "/var/foo"} # Missing key-name #2
    yield _do_test_atomiccmd__paths__malformed, {"TEMP_OUT"  : "/var/foo"} # Missing key-name #3
    yield _do_test_atomiccmd__paths__malformed, {"TEMP_OUT_" : "/var/foo"} # Missing key-name #4
    yield _do_test_atomiccmd__paths__malformed, {"TEMP_OUX_FOO" : "foo"} # Invalid key-type #1
    yield _do_test_atomiccmd__paths__malformed, {"INS_BAR"      : "foo"} # Invalid key-type #2


def test_atomiccmd__paths__invalid_values():
    @nose.tools.raises(TypeError)
    def _do_test_atomiccmd__paths__invalid_values(kwargs):
        AtomicCmd("true", **kwargs)

    yield _do_test_atomiccmd__paths__invalid_values, {"IN_FILE"       : 1}
    yield _do_test_atomiccmd__paths__invalid_values, {"TEMP_IN_FILE"  : set()}
    yield _do_test_atomiccmd__paths__invalid_values, {"OUT_FILE"      : [1, 2, 3]}
    yield _do_test_atomiccmd__paths__invalid_values, {"TEMP_OUT_FILE" : 1.0}

    yield _do_test_atomiccmd__paths__invalid_values, {"IN_STDIN"        : {}}
    yield _do_test_atomiccmd__paths__invalid_values, {"TEMP_IN_STDIN"   : frozenset()}
    yield _do_test_atomiccmd__paths__invalid_values, {"OUT_STDOUT"      : 1.7}
    yield _do_test_atomiccmd__paths__invalid_values, {"TEMP_OUT_STDOUT" : ()}
    yield _do_test_atomiccmd__paths__invalid_values, {"OUT_STDERR"      : xrange(3)}
    yield _do_test_atomiccmd__paths__invalid_values, {"TEMP_OUT_STDERR" : -1}


# Subpaths are not allowed for temp IN/OUT files, neither relative nor asbsolute
def test_atomiccmd__paths__invalid_temp_paths():
    @nose.tools.raises(ValueError)
    def _do_test_atomiccmd__paths__invalid_temp_paths(kwargs):
        AtomicCmd("true", **kwargs)

    # No relative paths
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_IN_FOO"     : "sub/infile"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_IN_STDIN"   : "sub/stdin"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_OUT_FOO"    : "sub/outfile"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_OUT_STDOUT" : "sub/stdout"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_OUT_STDERR" : "sub/stderr"}

    # No absolute paths
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_IN_FOO"     : "/tmp/sub/infile"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_IN_STDIN"   : "/dev/sub/stdin"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_OUT_FOO"    : "/etc/sub/outfile"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_OUT_STDOUT" : "/var/sub/stdout"}
    yield _do_test_atomiccmd__paths__invalid_temp_paths, {"TEMP_OUT_STDERR" : "/home/sub/stderr"}


# All OUT_ files must be unique, including all TEMP_OUT_
def test_atomiccmd__paths__overlapping_output():
    @nose.tools.raises(ValueError)
    def _do_test_atomiccmd__paths__overlapping_output(key_1, file_1, key_2, file_2):
        AtomicCmd(("ls",), **{key_1 : file_1, key_2 : file_2})

    yield _do_test_atomiccmd__paths__overlapping_output, "OUT_FILE_1", "/foo/bar/outfile", "OUT_FILE_2", "/var/outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "TEMP_OUT_FILE_1", "outfile", "OUT_FILE_1", "/var/outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "OUT_FILE_1", "/foo/bar/outfile", "TEMP_OUT_FILE_1", "outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "TEMP_OUT_FILE_1", "outfile", "TEMP_OUT_FILE_2", "outfile"

    yield _do_test_atomiccmd__paths__overlapping_output, "OUT_FILE_1", "/foo/bar/outfile", "OUT_STDOUT", "/var/outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "TEMP_OUT_FILE_1", "outfile", "OUT_STDOUT", "/var/outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "OUT_FILE_1", "/foo/bar/outfile", "TEMP_OUT_STDOUT", "outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "TEMP_OUT_FILE_1", "outfile", "TEMP_OUT_STDOUT", "outfile"

    yield _do_test_atomiccmd__paths__overlapping_output, "OUT_FILE_1", "/foo/bar/outfile", "OUT_STDERR", "/var/outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "TEMP_OUT_FILE_1", "outfile", "OUT_STDERR", "/var/outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "OUT_FILE_1", "/foo/bar/outfile", "TEMP_OUT_STDERR", "outfile"
    yield _do_test_atomiccmd__paths__overlapping_output, "TEMP_OUT_FILE_1", "outfile", "TEMP_OUT_STDERR", "outfile"


# A pipe can be w/wo TEMP_, but not both
def test_atomiccmd__pipes__duplicates():
    @nose.tools.raises(CmdError)
    def _do_test_atomiccmd__pipes__duplicates(key):
        kwargs = {"TEMP_" + key : "temp_file",
                  key           : "file"}
        AtomicCmd(["ls"], **kwargs)

    yield _do_test_atomiccmd__pipes__duplicates, "IN_STDIN"
    yield _do_test_atomiccmd__pipes__duplicates, "OUT_STDOUT"
    yield _do_test_atomiccmd__pipes__duplicates, "OUT_STDERR"



################################################################################
################################################################################
## CHECK_ / EXEC_

# RequirementObjs are the standard way to do tests
def test_atomicmcd__exec__reqobj():
    reqobj = RequirementObj(call = ("echo", "version"),
                            search = "version",
                            pprint = "{}",
                            checks = str)
    cmd = AtomicCmd("true",
                    CHECK_VERSION = reqobj)
    assert_equal(cmd.requirements, frozenset([reqobj]))

# CHECK_ is expected to be a callable
@nose.tools.raises(TypeError)
def test_atomiccmd__checks__non_callable():
    AtomicCmd("ls", CHECK_FOO = "ls")


# EXEC_ is expected to be a string
def test_atomiccmd__exec__invalid():
    @nose.tools.raises(TypeError)
    def _test_atomiccmd__exec__invalid(obj):
        AtomicCmd("true", EXEC_FOO = obj)

    yield _test_atomiccmd__exec__invalid, str
    yield _test_atomiccmd__exec__invalid, {}
    yield _test_atomiccmd__exec__invalid, 1


################################################################################
################################################################################
## Path components

@with_temp_folder
def test_atomiccmd__paths_non_str(temp_folder):
    cmd = AtomicCmd(("touch", 1234),
                    OUT_FOO = "1234",
                    set_cwd = True)
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    assert os.path.exists(os.path.join(temp_folder, "1234"))


@nose.tools.raises(CmdError)
def test_atomiccmd__paths_missing():
    AtomicCmd(("touch", "%(IN_FOO)s"))

@nose.tools.raises(CmdError)
def test_atomiccmd__paths_invalid():
    AtomicCmd(("touch", "%(IN_FOO)"),
               IN_FOO = "abc")


@with_temp_folder
def test_atomiccmd__paths__key(temp_folder):
    cmd = AtomicCmd(("echo", "-n", "%(TEMP_DIR)s"),
                    OUT_STDOUT = AtomicCmd.PIPE)
    cmd.run(temp_folder)
    path = cmd.stdout.read()
    assert os.path.samefile(temp_folder, path), (temp_folder, path)
    assert_equal(cmd.join(), [0])



################################################################################
################################################################################
## Constructor: Piping commands

@with_temp_folder
def test_atomiccmd__piping(temp_folder):
    cmd_1 = AtomicCmd(["echo", "-n", "#@!$^"],
                        OUT_STDOUT = AtomicCmd.PIPE)
    assert_equal(cmd_1.output_files, frozenset())
    cmd_2 = AtomicCmd(["cat"],
                        IN_STDIN = cmd_1,
                        OUT_STDOUT = "piped.txt")
    assert_equal(cmd_2.input_files, frozenset())
    cmd_1.run(temp_folder)
    cmd_2.run(temp_folder)
    assert_equal(cmd_1.join(), [0])
    assert_equal(cmd_2.join(), [0])
    result = get_file_contents(os.path.join(temp_folder, "piped.txt"))
    assert_equal(result, "#@!$^")


@with_temp_folder
def test_atomiccmd__piping_temp(temp_folder):
    cmd_1 = AtomicCmd(["echo", "-n", "#@!$^"],
                        TEMP_OUT_STDOUT = AtomicCmd.PIPE)
    assert_equal(cmd_1.output_files, frozenset())
    cmd_2 = AtomicCmd(["cat"],
                        TEMP_IN_STDIN = cmd_1,
                        OUT_STDOUT = "piped.txt")
    assert_equal(cmd_2.input_files, frozenset())
    cmd_1.run(temp_folder)
    cmd_2.run(temp_folder)
    assert_equal(cmd_1.join(), [0])
    assert_equal(cmd_2.join(), [0])
    result = get_file_contents(os.path.join(temp_folder, "piped.txt"))
    assert_equal(result, "#@!$^")

# Only STDOUT takes AtomicCmd.PIPE
def test_atomiccmd__piping__wrong_pipe():
    @nose.tools.raises(TypeError)
    def _test_atomiccmd__piping__wrong_pipe(key):
        AtomicCmd("ls", **{key : AtomicCmd.PIPE})

    yield _test_atomiccmd__piping__wrong_pipe, "IN_STDIN"
    yield _test_atomiccmd__piping__wrong_pipe, "TEMP_IN_STDIN"
    yield _test_atomiccmd__piping__wrong_pipe, "OUT_STDERR"
    yield _test_atomiccmd__piping__wrong_pipe, "TEMP_OUT_STDERR"




################################################################################
################################################################################
## run

@with_temp_folder
@nose.tools.raises(CmdError)
def test_atomiccmd__run__already_running(temp_files):
    cmd = AtomicCmd(("sleep", "10"))
    try:
        cmd.run(temp_files)
        cmd.run(temp_files)
    finally:
        cmd.terminate()



################################################################################
################################################################################
## Ready

@with_temp_folder
def test_atomiccmd__ready(temp_folder):
    cmd = AtomicCmd("ls")
    assert_equal(cmd.join(), [None])
    assert not cmd.ready()
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    assert cmd.ready()



################################################################################
################################################################################
## Join / wait

def test_atomiccmd__join_wait():
    @with_temp_folder
    def _do_test_atomiccmd__join_wait(temp_folder, func, call, before_run, after_run):
        cmd = AtomicCmd(call)
        assert_equal(func(cmd), before_run)
        cmd.run(temp_folder)
        assert_equal(func(cmd), after_run)

    yield _do_test_atomiccmd__join_wait, AtomicCmd.join, "true",  [None], [0]
    yield _do_test_atomiccmd__join_wait, AtomicCmd.join, "false", [None], [1]

    yield _do_test_atomiccmd__join_wait, AtomicCmd.wait, "true",  None,   0
    yield _do_test_atomiccmd__join_wait, AtomicCmd.wait, "false", None,   1



################################################################################
################################################################################
## Terminate

class _TerminateMock(AtomicCmd):
    class _ProcMock:
        def __init__(self, outer):
            self._outer = outer
        def terminate(self):
            self._outer.was_terminated = True
            if self._outer.raise_on_terminate:
                raise OSError("KABOOM!")
        def poll(self): # pylint: disable=R0201
            return None

    def __init__(self, raise_on_terminate):
        AtomicCmd.__init__(self, ("true"))
        self._proc = _TerminateMock._ProcMock(self)
        self.was_terminated = False
        self.raise_on_terminate = raise_on_terminate


def test_atomiccmd__terminate():
    def _do_test_atomiccmd__terminate(raise_on_terminate):
        cmdmock = _TerminateMock(raise_on_terminate)
        assert cmdmock._proc
        cmdmock.terminate()
        assert cmdmock.was_terminated
        # If the proc was finished, terminate does not clear it
        assert bool(cmdmock._proc) == raise_on_terminate
    yield _do_test_atomiccmd__terminate, False
    yield _do_test_atomiccmd__terminate, True


# Ensure that no OSException is raised, even if the command
# managed to finish before terminate was called
@with_temp_folder
def test_atomiccmd__terminate_race_condition(temp_folder):
    cmd = AtomicCmd("true")
    cmd.run(temp_folder)
    while cmd._proc.poll() is None:
        pass
    cmd.terminate()
    assert_equal(cmd.join(), [0])

# Calling terminate on an already joined command is acceptable ...
@with_temp_folder
def test_atomiccmd__terminate_after_join(temp_folder):
    cmd = AtomicCmd("true")
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    cmd.terminate()
    assert_equal(cmd.join(), [0])




################################################################################
################################################################################
## commit

def _setup_for_commit(temp_folder, create_cmd = True):
    destination = os.path.join(temp_folder, "out")
    temp_folder = os.path.join(temp_folder, "tmp")
    os.makedirs(destination)
    os.makedirs(temp_folder)

    if not create_cmd:
        return destination, temp_folder

    cmd = AtomicCmd(("touch", "%(OUT_FOO)s"),
                    OUT_FOO = os.path.join(destination, "1234"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])

    return destination, temp_folder, cmd


@with_temp_folder
def test_atomiccmd__commit_simple(temp_folder):
    destination, temp_folder, cmd = _setup_for_commit(temp_folder)
    cmd.commit(temp_folder)
    assert not os.path.exists(os.path.join(temp_folder, "1234"))
    assert os.path.exists(os.path.join(destination, "1234"))

@with_temp_folder
def test_atomiccmd__commit_temp_out(temp_folder):
    dest, temp = _setup_for_commit(temp_folder, create_cmd = False)
    cmd = AtomicCmd(("echo", "foo"),
                    OUT_STDOUT   = os.path.join(dest, "foo.txt"),
                    TEMP_OUT_FOO = "bar.txt")
    cmd.run(temp)
    assert_equal(cmd.join(), [0])
    set_file_contents(os.path.join(temp, "bar.txt"), "1 2 3")
    cmd.commit(temp)
    assert_equal(os.listdir(temp), [])
    assert_equal(os.listdir(dest), ["foo.txt"])

@with_temp_folder
def test_atomiccmd__commit_temp_only(temp_folder):
    cmd = AtomicCmd(("echo", "foo"),
                    TEMP_OUT_STDOUT = "bar.txt")
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    assert os.path.exists(os.path.join(temp_folder, "bar.txt"))
    cmd.commit(temp_folder)
    assert_equal(os.listdir(temp_folder), [])


@nose.tools.raises(CmdError)
def test_atomiccmd__commit_before_run():
    AtomicCmd("true").commit("/tmp")

@with_temp_folder
@nose.tools.raises(CmdError)
def test_atomiccmd__commit_while_running(temp_folder):
    cmd = AtomicCmd(("sleep", "10"))
    try:
        cmd.run(temp_folder)
        cmd.commit(temp_folder)
    finally:
        cmd.terminate()

@with_temp_folder
@nose.tools.raises(CmdError)
def test_atomiccmd__commit_before_join(temp_folder):
    cmd = AtomicCmd(("sleep", "0.1"))
    cmd.run(temp_folder)
    while cmd._proc.poll() is None:
        pass
    cmd.commit(temp_folder)

# The temp path might differ, as long as the actual path is the same
@with_temp_folder
def test_atomiccmd__commit_temp_folder(temp_folder):
    destination, temp_folder, cmd = _setup_for_commit(temp_folder)
    cmd.commit(os.path.realpath(temp_folder))
    assert not os.path.exists(os.path.join(temp_folder, "1234"))
    assert os.path.exists(os.path.join(destination, "1234"))


@with_temp_folder
@nose.tools.raises(CmdError)
def test_atomiccmd__commit_wrong_temp_folder(temp_folder):
    destination, temp_folder, cmd = _setup_for_commit(temp_folder)
    cmd.commit(destination)


@with_temp_folder
@nose.tools.raises(CmdError)
def test_atomiccmd__commit_missing_files(temp_folder):
    destination, temp_folder = _setup_for_commit(temp_folder, False)
    cmd = AtomicCmd(("touch", "%(OUT_FOO)s"),
                    OUT_FOO = os.path.join(destination, "1234"),
                    OUT_BAR = os.path.join(destination, "4567"))
    cmd.run(temp_folder)
    cmd.join()
    before = set(os.listdir(temp_folder))
    try:
        cmd.commit(temp_folder)
    except CmdError: # No changes should have been made yet
        assert before == set(os.listdir(temp_folder))
        raise




################################################################################
################################################################################
## stdout

def test_atomiccmd__stdout():
    @with_temp_folder
    def _do_test_atomiccmd__stdout(temp_folder, kwargs):
        cmd = AtomicCmd(("echo", "-n", "foo"), **kwargs)
        assert cmd.stdout is None
        cmd.run(temp_folder)
        assert isinstance(cmd.stdout, types.FileType), cmd.stdout
        assert_equal(cmd.stdout.read(), "foo")
        cmd.join()
        cmd.commit(temp_folder)
        assert cmd.stdout is None

    yield _do_test_atomiccmd__stdout, {"OUT_STDOUT" : AtomicCmd.PIPE}
    yield _do_test_atomiccmd__stdout, {"TEMP_OUT_STDOUT" : AtomicCmd.PIPE}

def test_atomiccmd__stdout_no_pipe():
    @with_temp_folder
    def _do_test_atomiccmd__stdout_no_pipe(temp_folder, kwargs):
        cmd = AtomicCmd(("echo", "foo"), **kwargs)
        assert cmd.stdout is None
        cmd.run(temp_folder)
        assert cmd.stdout is None
        cmd.join()

    yield _do_test_atomiccmd__stdout_no_pipe, {"OUT_STDOUT" : "foo.txt"}
    yield _do_test_atomiccmd__stdout_no_pipe, {"TEMP_OUT_STDOUT" : "foo.txt"}




################################################################################
################################################################################
## __str__
## Additional tests in atomicpp_test.py

def test_atomiccmd__str__():
    cmd = AtomicCmd(("echo", "test"))
    assert_equal(pypeline.atomicpp.pformat(cmd), str(cmd))



################################################################################
################################################################################
## Cleanup

# Test that the internal list of processes is kept clean of old objects
def test_atomiccmd__cleanup_proc():
    @with_temp_folder
    def _do_test_atomiccmd__cleanup_proc(temp_folder, func):
        assert_equal(pypeline.atomiccmd._PROCS, set())
        cmd = AtomicCmd("ls")
        cmd.run(temp_folder)
        ref = iter(pypeline.atomiccmd._PROCS).next()
        assert ref
        assert_equal(ref(), cmd._proc)

        assert_equal(cmd.join(), [0])
        cmd = func(cmd, temp_folder)

        assert ref not in pypeline.atomiccmd._PROCS

    def _do_commit(cmd, temp_folder):
        # Trigger freeing of proc
        cmd.commit(temp_folder)
        return cmd

    # The proc object should be released when commit is called
    yield _do_test_atomiccmd__cleanup_proc, _do_commit
    # The proc object should be released when the cmd object is released
    yield _do_test_atomiccmd__cleanup_proc, lambda _cmd, _temp_folder: None


# Cleanup of process object should be done if terminate is run
@nose.tools.timed(0.1)
@with_temp_folder
def test_atomiccmd__cleanup_proc__terminate(temp_folder):
    assert_equal(pypeline.atomiccmd._PROCS, set())
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(temp_folder)
    ref = iter(pypeline.atomiccmd._PROCS).next()
    assert ref
    assert_equal(ref(), cmd._proc)
    cmd.terminate()
    assert not cmd._proc
    assert not ref()
    assert ref not in pypeline.atomiccmd._PROCS


@with_temp_folder
def test_atomiccmd__cleanup_sigterm(temp_folder):
    sigs_sent, exit_called = {}, []
    def _wrap_killpg(pid, sig):
        assert pid not in sigs_sent
        sigs_sent[pid] = sig
    def _wrap_exit(rc):
        exit_called.append(rc)

    assert not pypeline.atomiccmd._PROCS
    cmd_1 = AtomicCmd(("sleep", "10"))
    cmd_2 = AtomicCmd(("sleep", "10"))
    cmd_1.run(temp_folder)
    cmd_2.run(temp_folder)

    with monkeypatch("os.killpg", _wrap_killpg):
        with monkeypatch("sys.exit", _wrap_exit):
            pypeline.atomiccmd._cleanup_children(signal.SIGTERM, None)

    assert_equal(exit_called, [-signal.SIGTERM])
    assert_equal(len(pypeline.atomiccmd._PROCS), 2)
    assert_equal(sigs_sent.get(cmd_1._proc.pid), signal.SIGTERM)
    assert_equal(sigs_sent.get(cmd_2._proc.pid), signal.SIGTERM)


# Ensure that the cleanup function handles weakrefs that have been freed
def test_atomiccmd__cleanup_sigterm__dead_weakrefs():
    exit_called = []
    procs_wrapper = [weakref.ref(monkeypatch("sys.exit", None))]

    assert_equal(procs_wrapper[0](), None)
    def _wrap_killpg(_pid, _sig):
        assert False # pragma: no coverage
    def _wrap_exit(rc):
        exit_called.append(rc)

    with monkeypatch("pypeline.atomiccmd._PROCS", procs_wrapper):
        with monkeypatch("os.killpg", _wrap_killpg):
            with monkeypatch("sys.exit", _wrap_exit):
                pypeline.atomiccmd._cleanup_children(signal.SIGTERM, None)
    assert_equal(exit_called, [-signal.SIGTERM])

