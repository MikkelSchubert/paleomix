#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is herby granted, free of charge, to any person obtaining a copy
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
# Disable warning for missign docstrings
# pylint: disable=C0111
# Disable warning caused by "invalid" function names
# pylint: disable=C0103
# Disable warning caused by touching private member variables/functions
# pylint: disable=W0212
# Disable warnings caused by flexmock setups ("X is assigned to nothing")
# pylint: disable=W0106
import os
import random

from nose.tools import \
    assert_in, \
    assert_equal, \
    assert_raises
from flexmock import flexmock

from paleomix.common.testing import \
    with_temp_folder, \
    set_file_contents, \
    get_file_contents

from paleomix.atomiccmd.command import \
    AtomicCmd
from paleomix.node import \
    Node, \
    CommandNode, \
    NodeError, \
    NodeUnhandledException, \
    CmdNodeError
from paleomix.common.utilities import \
    safe_coerce_to_frozenset


def test_dir():
    return os.path.dirname(__file__)


def test_file(*args):
    return os.path.join(test_dir(), "data", *args)


def _CommandNodeWrap(**kwargs):
    return CommandNode(command=AtomicCmd("true"), **kwargs)
_NODE_TYPES = (Node, _CommandNodeWrap)


_DESCRIPTION = "My description of a node"
_IN_FILES = frozenset((test_file("empty_file_1"),
                       test_file("empty_file_2")))
_OUT_FILES = frozenset((test_file("missing_out_file_1"),
                        test_file("missing_out_file_2")))
_EXEC_FILES = frozenset(("ls", "sh"))
_AUX_FILES = frozenset((test_file("rCRS.fasta"),
                        test_file("rCRS.fasta.fai")))
_REQUIREMENTS = frozenset((id, str))
_EMPTY_FILE = test_file("empty_file_1")


def _build_cmd_mock(input_files  = _IN_FILES,
                    output_files = (),
                    executables  = (),
                    auxiliary_files = (),
                    requirements = (),
                    optional_temp_files = ()):
    return flexmock(input_files     = frozenset(input_files),
                    output_files    = frozenset(output_files),
                    executables     = frozenset(executables),
                    auxiliary_files = frozenset(auxiliary_files),
                    requirements    = frozenset(requirements),
                    expected_temp_files = frozenset(map(os.path.basename, output_files)),
                    optional_temp_files = frozenset(optional_temp_files))


###############################################################################
###############################################################################
# Node: Constructor: File sets

def test_constructor():
    def first(values):
        return random.choice(tuple(values))

    def _do_test_constructor__single_value(key, value):
        defaults = {"input_files": _EMPTY_FILE}
        defaults[key] = value
        node = Node(**defaults)
        expected = safe_coerce_to_frozenset(value)
        assert_equal(getattr(node, key), expected)

    # Single values
    yield _do_test_constructor__single_value, "input_files",  first(_IN_FILES)
    yield _do_test_constructor__single_value, "output_files", first(_OUT_FILES)
    yield _do_test_constructor__single_value, "executables",  first(_EXEC_FILES)
    yield _do_test_constructor__single_value, "auxiliary_files",  first(_AUX_FILES)

    # Single value in list
    yield _do_test_constructor__single_value, "input_files",  [first(_IN_FILES)]
    yield _do_test_constructor__single_value, "output_files", [first(_OUT_FILES)]
    yield _do_test_constructor__single_value, "executables",  [first(_EXEC_FILES)]
    yield _do_test_constructor__single_value, "auxiliary_files",  [first(_AUX_FILES)]

    # Multiple values in list
    yield _do_test_constructor__single_value, "input_files",  _IN_FILES
    yield _do_test_constructor__single_value, "output_files", _OUT_FILES
    yield _do_test_constructor__single_value, "executables",  _EXEC_FILES
    yield _do_test_constructor__single_value, "auxiliary_files",  _AUX_FILES


def test_constructor__invalid_values():
    def _do_test_constructor__invalid_values(key, value):
        assert_raises(TypeError, Node, **{key: value})

    yield _do_test_constructor__invalid_values, "input_files",      [id]
    yield _do_test_constructor__invalid_values, "output_files",     [-1]
    yield _do_test_constructor__invalid_values, "executables",      [{}]
    yield _do_test_constructor__invalid_values, "auxiliary_files",  [1.3]


###############################################################################
###############################################################################
# Node: Constructor: Requirements

def test_constructor__requirements():
    node = Node(requirements=id)
    assert_equal(node.requirements, frozenset([id]))
    node = Node(requirements=[id])
    assert_equal(node.requirements, frozenset([id]))
    node = Node(requirements=[id, str])
    assert_equal(node.requirements, frozenset([id, str]))


def test_constructor__requirements__wrong_type():
    def _do_test_constructor__requirements__wrong_type(value):
        assert_raises(TypeError, Node, requirements=value)

    yield _do_test_constructor__requirements__wrong_type, 17
    yield _do_test_constructor__requirements__wrong_type, {}
    yield _do_test_constructor__requirements__wrong_type, "867-5309"


###############################################################################
###############################################################################
# Node: Constructor: Dependencies

def test_constructor__nodes_is_none():
    my_node = Node(dependencies=None)
    assert_equal(my_node.dependencies, frozenset())


def test_constructor__single_node():
    sub_node = Node()
    my_node = Node(dependencies=sub_node)
    assert_equal(my_node.dependencies, frozenset([sub_node]))


def test_constructor__iterable():
    sub_nodes = [Node(), Node()]
    my_node = Node(dependencies=iter(sub_nodes))
    assert_equal(my_node.dependencies, frozenset(sub_nodes))


def test_constructor__not_a_node():
    assert_raises(TypeError, Node, dependencies=(1,))


###############################################################################
###############################################################################
# *Node: Description

def test_constructor__description():
    def _do_test_constructor__description(cls):
        my_node = cls(description=_DESCRIPTION)
        assert_equal(str(my_node), _DESCRIPTION)
    for cls in _NODE_TYPES:
        yield _do_test_constructor__description, cls


def test_constructor__description__default():
    def _do_test_constructor__description__default(cls):
        my_node = cls()
        assert_equal(str(my_node), repr(my_node))
    for cls in _NODE_TYPES:
        yield _do_test_constructor__description__default, cls


def test_constructor__description__non_string():
    def _do_test_constructor__description__non_string(cls, value):
        assert_raises(TypeError, cls, description=value)
    for cls in _NODE_TYPES:
        yield _do_test_constructor__description__non_string, cls, 1
        yield _do_test_constructor__description__non_string, cls, {}


###############################################################################
###############################################################################
# *Node: Constructor tests: #threads

def test_constructor__threads():
    def _do_test_constructor__threads(cls, nthreads):
        node = cls(threads=nthreads)
        assert_equal(node.threads, nthreads)
    for cls in (Node, _CommandNodeWrap):
        yield _do_test_constructor__threads, cls, 1L
        yield _do_test_constructor__threads, cls, 3


def test_constructor__threads_invalid_range():
    def _do_test_constructor__threads_invalid_range(cls, nthreads):
        assert_raises(ValueError, cls, threads=nthreads)
    for cls in (Node, _CommandNodeWrap):
        yield _do_test_constructor__threads_invalid_range, cls, -1
        yield _do_test_constructor__threads_invalid_range, cls, 0


def test_constructor__threads_invalid_type():
    def _do_test_constructor__threads_invalid_type(cls, nthreads):
        assert_raises(TypeError, cls, threads=nthreads)
    for cls in (Node, _CommandNodeWrap):
        yield _do_test_constructor__threads_invalid_type, cls, "1"
        yield _do_test_constructor__threads_invalid_type, cls, {}
        yield _do_test_constructor__threads_invalid_type, cls, 2.7


###############################################################################
###############################################################################
# Node: Run

_DUMMY_TEMP_ROOT = "/xyz/tmp"
_DUMMY_TEMP = os.path.join(_DUMMY_TEMP_ROOT, "xTMPx")


def test_run__order():
    cfg_mock = flexmock(temp_root=_DUMMY_TEMP_ROOT)
    node_mock = flexmock(Node())
    node_mock.should_receive("_create_temp_dir").with_args(cfg_mock).and_return(_DUMMY_TEMP).ordered.once
    node_mock.should_receive("_setup").with_args(cfg_mock, _DUMMY_TEMP).ordered.once
    node_mock.should_receive("_run").with_args(cfg_mock, _DUMMY_TEMP).ordered.once
    node_mock.should_receive("_teardown").with_args(cfg_mock, _DUMMY_TEMP).ordered.once
    node_mock.should_receive("_remove_temp_dir").with_args(_DUMMY_TEMP).ordered.once
    node_mock.run(cfg_mock)


def test_run__exceptions():
    cfg_mock = flexmock(temp_root=_DUMMY_TEMP_ROOT)

    def build_tests(key, exception, expectation):
        def test_function():
            node_mock = flexmock(Node())
            node_mock.should_receive('_create_temp_dir').with_args(cfg_mock) \
                .and_return(_DUMMY_TEMP).ordered.once
            node_mock.should_receive(key).and_raise(exception).ordered.once
            node_mock.should_receive('_remove_temp_dir').never

            assert_raises(expectation, node_mock.run, cfg_mock)

        return test_function

    print "foo"
    for key in ('_setup', '_run', '_teardown'):
        yield build_tests(key, TypeError("The castle AAARGH!"), NodeUnhandledException)
        yield build_tests(key, NodeError("He's a very naughty boy!"), NodeError)


def test_run__exception__create_temp_dir():
    cfg_mock = flexmock(temp_root=_DUMMY_TEMP_ROOT)
    node_mock = flexmock(Node())
    node_mock.should_receive('_create_temp_dir').with_args(cfg_mock) \
        .and_raise(OSError()).ordered.once

    assert_raises(NodeUnhandledException, node_mock.run, cfg_mock)


def test_run__exception__remove_temp_dir():
    cfg_mock = flexmock(temp_root=_DUMMY_TEMP_ROOT)
    node_mock = flexmock(Node())
    node_mock.should_receive('_create_temp_dir').with_args(cfg_mock) \
        .and_return(_DUMMY_TEMP).ordered.once
    node_mock.should_receive('_remove_temp_dir').with_args(_DUMMY_TEMP) \
        .and_raise(OSError()).ordered.once

    assert_raises(NodeUnhandledException, node_mock.run, cfg_mock)


def test_run__error_log__node_error():
    @with_temp_folder
    def _do_test_run__error_log__node_error(temp_folder, exception):
        temp = os.path.join(temp_folder, "xTMPx")
        cfg_mock = flexmock(temp_root=temp_folder)
        node_mock = flexmock(Node())
        node_mock.should_receive("_create_temp_dir").with_args(cfg_mock) \
            .and_return(temp).ordered.once
        node_mock.should_receive("_run").and_raise(exception).ordered.once

        os.mkdir(temp)
        assert_raises(NodeError, node_mock.run, cfg_mock)
        log_file = os.path.join(temp_folder, "xTMPx", "pipe.errors")
        assert os.path.exists(log_file)
        assert_in("Errors =", get_file_contents(log_file))

    yield _do_test_run__error_log__node_error, NodeError("ARGH!")
    yield _do_test_run__error_log__node_error, OSError("ARGH!")


###############################################################################
###############################################################################
# Node: _setup / _teardown

def test__setup__input_files():
    def _do_test__setup__input_files_exist(kwargs):
        Node(**kwargs)._setup(None, None)

    yield _do_test__setup__input_files_exist, {"executables": ("ls", "sh")}
    yield _do_test__setup__input_files_exist, {"input_files": _IN_FILES}
    yield _do_test__setup__input_files_exist, {"auxiliary_files": _IN_FILES}


def test__setup__input_files_missing():
    def _do_test__setup__input_files_exist(kwargs):
        assert_raises(NodeError, Node(**kwargs)._setup, None, None)

    yield _do_test__setup__input_files_exist, {"executables": ("ls", "shxxxx")}
    yield _do_test__setup__input_files_exist, {"input_files": _OUT_FILES}
    yield _do_test__setup__input_files_exist, {"auxiliary_files": _OUT_FILES}


def test__teardown__output_files():
    Node(input_files=_EMPTY_FILE,
         output_files=_IN_FILES)._teardown(None, None)


def test__teardown__output_files_missing():
    node = Node(input_files=_EMPTY_FILE,
                output_files=_OUT_FILES)
    assert_raises(NodeError, node._teardown, None, None)


###############################################################################
###############################################################################
# CommandNode: Constructor

_SIMPLE_DEPS = Node()
_SIMPLE_SUBS = Node()
_SIMPLE_CMD_MOCK = flexmock(input_files=_IN_FILES,
                            output_files=_OUT_FILES,
                            executables=_EXEC_FILES,
                            auxiliary_files=_AUX_FILES,
                            requirements=_REQUIREMENTS)
_SIMPLE_CMD_NODE = CommandNode(command=_SIMPLE_CMD_MOCK,
                               dependencies=_SIMPLE_DEPS)


def test_commandnode_constructor__input_files():
    assert_equal(_SIMPLE_CMD_NODE.input_files, _IN_FILES)


def test_commandnode_constructor__output_files():
    assert_equal(_SIMPLE_CMD_NODE.output_files, _OUT_FILES)


def test_commandnode_constructor__auxiliary_files():
    assert_equal(_SIMPLE_CMD_NODE.auxiliary_files, _AUX_FILES)


def test_commandnode_constructor__executables():
    assert_equal(_SIMPLE_CMD_NODE.executables, _EXEC_FILES)


def test_commandnode_constructor__requirements():
    assert_equal(_SIMPLE_CMD_NODE.requirements, _REQUIREMENTS)


def test_commandnode_constructor__dependencies():
    assert_equal(_SIMPLE_CMD_NODE.dependencies, frozenset([_SIMPLE_DEPS]))


def test_commandnode_constructor__dependencies__default():
    cmd_mock = CommandNode(command=_SIMPLE_CMD_MOCK)
    assert_equal(cmd_mock.dependencies, frozenset())


###############################################################################
###############################################################################
# CommandNode: run

def test_command_node__run():
    cfg_mock = flexmock(temp_root=_DUMMY_TEMP_ROOT)
    cmd_mock = _build_cmd_mock()
    node_mock = flexmock(CommandNode(cmd_mock))
    node_mock.should_receive("_create_temp_dir").with_args(cfg_mock) \
        .and_return(_DUMMY_TEMP).ordered.once
    node_mock.should_receive("_setup").with_args(cfg_mock, _DUMMY_TEMP).ordered.once
    cmd_mock.should_receive("run").with_args(_DUMMY_TEMP).ordered.once
    cmd_mock.should_receive("join").and_return([0]).ordered.once
    node_mock.should_receive("_teardown").with_args(cfg_mock, _DUMMY_TEMP).ordered.once
    node_mock.should_receive("_remove_temp_dir").with_args(_DUMMY_TEMP).ordered.once
    node_mock.run(cfg_mock)


###############################################################################
###############################################################################
# CommandNode: _setup()

def test_commandnode_setup__files_exist():
    def _do_test_commandnode_setup(kwargs):
        cmd_mock = _build_cmd_mock(**kwargs)
        node = CommandNode(cmd_mock)
        node._setup(None, None)
    yield _do_test_commandnode_setup, {"executables": ("ls", "sh")}
    yield _do_test_commandnode_setup, {"input_files": _IN_FILES}
    yield _do_test_commandnode_setup, {"auxiliary_files": _IN_FILES}


def test_commandnode_setup__files_missing():
    def _do_test_commandnode_setup(kwargs):
        cmd_mock = _build_cmd_mock(**kwargs)
        node = CommandNode(cmd_mock)
        assert_raises(NodeError, node._setup, None, None)

    yield _do_test_commandnode_setup, {"executables": ("ls", "shxxxxxxxxxxx")}
    yield _do_test_commandnode_setup, {"input_files": _OUT_FILES}
    yield _do_test_commandnode_setup, {"auxiliary_files": _OUT_FILES}


###############################################################################
###############################################################################
# CommandNode: _run()

def test_commandnode_run__call_order():
    cmd_mock = _build_cmd_mock()
    cmd_mock.should_receive("run").with_args("xTMPx").ordered.once
    cmd_mock.should_receive("join").with_args().and_return((0,)).ordered.once
    node = CommandNode(cmd_mock)
    node._run(None, "xTMPx")


def test_commandnode_run__exception_on_error():
    cmd_mock = _build_cmd_mock()
    cmd_mock.should_receive("run").ordered.once
    cmd_mock.should_receive("join").and_return((1,)).ordered.once
    node = CommandNode(cmd_mock)
    assert_raises(CmdNodeError, node._run, None, None)


###############################################################################
###############################################################################
# CommandNode: _teardown

def _setup_temp_folders(temp_folder):
    destination = os.path.join(temp_folder, "dst")
    temp_folder = os.path.join(temp_folder, "tmp")
    os.makedirs(temp_folder)
    os.makedirs(destination)
    return destination, temp_folder


# Commit is called on the command obj
@with_temp_folder
def test_commandnode_teardown__commit(temp_folder):
    cmd_mock = _build_cmd_mock()
    cmd_mock.should_receive("commit").with_args(temp_folder).once
    node = CommandNode(cmd_mock)
    node._teardown(None, temp_folder)


# Files exist in temp folder, and in destination after commit
@with_temp_folder
def test_commandnode_teardown(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    cmd = AtomicCmd(("echo", "-n", "1 2 3"),
                    IN_DUMMY=_EMPTY_FILE,
                    OUT_STDOUT=os.path.join(destination, "foo.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    assert os.path.exists(os.path.join(temp_folder, "foo.txt"))
    assert not os.path.exists(os.path.join(destination, "foo.txt"))
    node._teardown(None, temp_folder)
    assert not os.path.exists(os.path.join(temp_folder, "foo.txt"))
    assert os.path.exists(os.path.join(destination, "foo.txt"))


# Not all required files have been generated (atomic)
@with_temp_folder
def test_commandnode_teardown__missing_files_in_temp(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    cmd = AtomicCmd(("echo", "-n", "1 2 3"),
                    IN_DUMMY=_EMPTY_FILE,
                    OUT_BAR=os.path.join(destination, "bar.txt"),
                    OUT_STDOUT=os.path.join(destination, "foo.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    temp_files_before = set(os.listdir(temp_folder))
    dest_files_before = set(os.listdir(destination))

    assert_raises(CmdNodeError, node._teardown, None, temp_folder)
    assert_equal(temp_files_before, set(os.listdir(temp_folder)))
    assert_equal(dest_files_before, set(os.listdir(destination)))


# Not all specified TEMP_ files exist at _teardown (allowed)
@with_temp_folder
def test_commandnode_teardown__missing_optional_files(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    cmd = AtomicCmd(("echo", "-n", "1 2 3"),
                    IN_DUMMY=_EMPTY_FILE,
                    TEMP_OUT_BAR="bar.txt",
                    OUT_STDOUT=os.path.join(destination, "foo.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    node._teardown(None, temp_folder)
    assert_equal(os.listdir(temp_folder), [])
    assert_equal(os.listdir(destination), ["foo.txt"])


# Not all required files were in place after commit
@with_temp_folder
def _test_commandnode_teardown__missing_files_in_dest(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    class _CmdMock(AtomicCmd):
        def commit(self, temp):
            AtomicCmd.commit(self, temp)
            os.remove(os.path.join(destination, "foo.txt"))

    cmd = _CmdMock(("touch", "%(OUT_FOO)s", "%(OUT_BAR)s"),
                   IN_DUMMY=_EMPTY_FILE,
                   OUT_FOO=os.path.join(destination, "foo.txt"),
                   OUT_BAR=os.path.join(destination, "bar.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    assert_raises(NodeError, node._teardown, None, temp_folder)


# Unexpected files were found in the temporary directory
@with_temp_folder
def test_commandnode_teardown__extra_files_in_temp(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    cmd = AtomicCmd(("echo", "-n", "1 2 3"),
                    IN_DUMMY=_EMPTY_FILE,
                    OUT_STDOUT=os.path.join(destination, "foo.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    set_file_contents(os.path.join(temp_folder, "bar.txt"), "1 2 3")
    temp_files_before = set(os.listdir(temp_folder))
    dest_files_before = set(os.listdir(destination))

    assert_raises(CmdNodeError, node._teardown, None, temp_folder)
    assert_equal(temp_files_before, set(os.listdir(temp_folder)))
    assert_equal(dest_files_before, set(os.listdir(destination)))
