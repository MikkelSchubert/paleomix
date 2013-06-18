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
# Disable warning caused by "invalid" function names
# pylint: disable=C0103
# Disable warning caused by touching private member variables/functions
# pylint: disable=W0212
# Disable warnings caused by flexmock setups ("X is assigned to nothing")
# pylint: disable=W0106
import os

import nose.tools
from nose.tools import assert_equal, assert_in # pylint: disable=E0611
from flexmock import flexmock
from tests.common.utils import monkeypatch, with_temp_folder, \
     set_file_contents, \
     get_file_contents

from pypeline.atomiccmd import AtomicCmd
from pypeline.node import Node, MetaNode, CommandNode, NodeError, NodeUnhandledException, \
     MetaNodeError, CmdNodeError
from pypeline.common.utilities import safe_coerce_to_frozenset




def _CommandNodeWrap(**kwargs):
    return CommandNode(command = AtomicCmd("true"), **kwargs)
_NODE_TYPES = (Node, _CommandNodeWrap, MetaNode)


class monkeypatch_tmp_dir:
    """Monkeypatches functions used by Node.run to setup the
    temporary folders. This is done to reduce the amount of
    file operations that actually need to be done."""
    def __init__(self, root = "/tmp", subfolder = "xTMPx"):
        self._mkdir_patch = monkeypatch("pypeline.common.fileutils.create_temp_dir", self._mkdir)
        self._rmdir_patch = monkeypatch("os.rmdir", self._rmdir)
        self._root_dir = root
        self._sub_dir  = subfolder
        self._mkdir_called = False
        self._rmdir_called = False

    def __enter__(self):
        self._mkdir_patch.__enter__()
        self._rmdir_patch.__enter__()
        return self

    def __exit__(self, type, value, traceback): # pylint: disable=W0622
        self._mkdir_patch.__exit__(type, value, traceback)
        self._rmdir_patch.__exit__(type, value, traceback)
        if not value:
            assert self._mkdir_called
            assert self._rmdir_called

    def _mkdir(self, path):
        self._mkdir_called = True
        assert_equal(self._root_dir, path)
        return os.path.join(self._root_dir, self._sub_dir)

    def _rmdir(self, path):
        self._rmdir_called = True
        return assert_equal(os.path.join(self._root_dir, self._sub_dir), path)



_DESCRIPTION = "My description of a node"
_IN_FILES    = frozenset(("tests/data/empty_file_1", "tests/data/empty_file_2"))
_OUT_FILES   = frozenset(("tests/data/missing_out_file_1", "tests/data/missing_out_file_2"))
_EXEC_FILES  = frozenset(("ls", "sh"))
_AUX_FILES   = frozenset(("tests/data/rCRS.fasta", "tests/data/rCRS.fasta.fai"))
_REQUIREMENTS = frozenset((id, str))


def _build_cmd_mock(input_files  = (), output_files = (), executables  = (), auxiliary_files = (), requirements = (),
                    optional_temp_files = ()):
    return flexmock(input_files     = frozenset(input_files),
                    output_files    = frozenset(output_files),
                    executables     = frozenset(executables),
                    auxiliary_files = frozenset(auxiliary_files),
                    requirements    = frozenset(requirements),
                    expected_temp_files = frozenset(map(os.path.basename, output_files)),
                    optional_temp_files = frozenset(optional_temp_files))



################################################################################
################################################################################
## Node: Constructor: File sets

def test_constructor():
    import random
    def first(values):
        return random.choice(tuple(values))

    def _do_test_constructor__single_value(key, value):
        node   = Node(**{key : value})
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
    @nose.tools.raises(TypeError)
    def _do_test_constructor__invalid_values(key, value):
        Node(**{key : value})

    yield _do_test_constructor__invalid_values, "input_files",      [id]
    yield _do_test_constructor__invalid_values, "output_files",     [-1]
    yield _do_test_constructor__invalid_values, "executables",      [{}]
    yield _do_test_constructor__invalid_values, "auxiliary_files",  [1.3]



################################################################################
################################################################################
## Node: Constructor: Requirements

def test_constructor__requirements():
    node = Node(requirements = id)
    assert_equal(node.requirements, frozenset([id]))
    node = Node(requirements = [id])
    assert_equal(node.requirements, frozenset([id]))
    node = Node(requirements = [id, str])
    assert_equal(node.requirements, frozenset([id, str]))

def test_constructor__requirements__wrong_type():
    @nose.tools.raises(TypeError)
    def _do_test_constructor__requirements__wrong_type(value):
        Node(requirements = value)
    yield _do_test_constructor__requirements__wrong_type, 17
    yield _do_test_constructor__requirements__wrong_type, {}
    yield _do_test_constructor__requirements__wrong_type, "867-5309"

@nose.tools.raises(NodeError)
def test_constructor__requirements__non_picklable():
    unpicklable = lambda: None # pragma: no coverage
    Node(requirements = unpicklable)


class _UnpicklableNode(Node):
    def __init__(self):
        self.a_property = lambda: None # pragma: no coverage
        Node.__init__(self)

@nose.tools.raises(NodeError)
def test_constructor__requirements__non_picklable_property():
    _UnpicklableNode()



################################################################################
################################################################################
## Node: Constructor: Dependencies / Subnodes

def test_constructor__nodes_is_none():
    def _do_test_constructor__none_nodes(key):
        my_node = Node(**{key : None})
        assert_equal(getattr(my_node, key), frozenset())
    yield _do_test_constructor__none_nodes, "subnodes"
    yield _do_test_constructor__none_nodes, "dependencies"

def test_constructor__single_node():
    def _do_test_constructor__single_node(key):
        sub_node = Node()
        my_node  = Node(**{key : sub_node})
        assert_equal(getattr(my_node, key), frozenset([sub_node]))
    yield _do_test_constructor__single_node, "subnodes"
    yield _do_test_constructor__single_node, "dependencies"

def test_constructor__iterable():
    def _do_test_constructor__iterable(key):
        sub_nodes = [Node(), Node()]
        my_node  = Node(**{key : iter(sub_nodes)})
        assert_equal(getattr(my_node, key), frozenset(sub_nodes))
    yield _do_test_constructor__iterable, "subnodes"
    yield _do_test_constructor__iterable, "dependencies"

def test_constructor__not_a_node():
    @nose.tools.raises(TypeError)
    def _do_test_constructor__not_a_node(key):
        Node(**{key : (1,)})
    yield _do_test_constructor__not_a_node, "subnodes"
    yield _do_test_constructor__not_a_node, "dependencies"



################################################################################
################################################################################
## *Node: Description


def test_constructor__description():
    def _do_test_constructor__description(cls):
        my_node = cls(description  = _DESCRIPTION)
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
    @nose.tools.raises(TypeError)
    def _do_test_constructor__description__non_string(cls, value):
        cls(description = value)
    for cls in _NODE_TYPES:
        yield _do_test_constructor__description__non_string, cls, 1
        yield _do_test_constructor__description__non_string, cls, {}




################################################################################
################################################################################
## *Node: Constructor tests: #threads

def test_constructor__threads():
    def _do_test_constructor__threads(cls, nthreads):
        node = cls(threads = nthreads)
        assert_equal(node.threads, nthreads)
    for cls in (Node, _CommandNodeWrap):
        yield _do_test_constructor__threads, cls, 1L
        yield _do_test_constructor__threads, cls, 3

def test_constructor__threads_invalid_range():
    @nose.tools.raises(ValueError)
    def _do_test_constructor__threads_invalid_range(cls, nthreads):
        cls(threads = nthreads)
    for cls in (Node, _CommandNodeWrap):
        yield _do_test_constructor__threads_invalid_range, cls, -1
        yield _do_test_constructor__threads_invalid_range, cls, 0

def test_constructor__threads_invalid_type():
    @nose.tools.raises(TypeError)
    def _do_test_constructor__threads_invalid_type(cls, nthreads):
        cls(threads = nthreads)
    for cls in (Node, _CommandNodeWrap):
        yield _do_test_constructor__threads_invalid_type, cls, "1"
        yield _do_test_constructor__threads_invalid_type, cls, {}
        yield _do_test_constructor__threads_invalid_type, cls, 2.7

def test_constuctor__threads_meta():
    node = MetaNode()
    assert_equal(node.threads, 1)




################################################################################
################################################################################
## Node: States

def test_is_done__no_output():
    assert Node().is_done

@with_temp_folder
def test_is_done__output_changes(temp_folder):
    temp_file_1 = os.path.join(temp_folder, "file_1.txt")
    temp_file_2 = os.path.join(temp_folder, "file_2.txt")
    my_node   = Node(output_files = (temp_file_1, temp_file_2))
    assert not my_node.is_done
    set_file_contents(temp_file_1, "foo")
    assert not my_node.is_done
    set_file_contents(temp_file_2, "bar")
    assert my_node.is_done

@with_temp_folder
def test_is_done__subnode_output_changes(temp_folder):
    temp_file = os.path.join(temp_folder, "file.txt")
    subnode   = Node(output_files = temp_file)
    my_node   = Node(subnodes = subnode)
    assert not my_node.is_done
    set_file_contents(temp_file, "foo")
    assert my_node.is_done


def test_is_outdated__no_output():
    assert not Node().is_outdated

def test_is_outdated__input_but_no_output():
    assert not Node(input_files = _IN_FILES).is_outdated

def test_is_outdated__output_but_no_input():
    assert not Node(output_files = _OUT_FILES).is_outdated

def test_is_outdated__not_outdated():
    my_node = Node(input_files  = "tests/data/timestamp_a_older",
                   output_files = "tests/data/timestamp_a_younger")
    assert not my_node.is_outdated

def test_is_outdated__outdated():
    my_node = Node(input_files  = "tests/data/timestamp_a_younger",
                   output_files = "tests/data/timestamp_a_older")
    assert my_node.is_outdated

def test_is_outdated__updates():
    my_node = Node(input_files  = "tests/data/timestamp_a_older",
                   output_files = "tests/data/timestamp_a_younger")
    assert not my_node.is_outdated
    my_node = Node(input_files  = "tests/data/timestamp_a_younger",
                   output_files = "tests/data/timestamp_a_older")
    assert my_node.is_outdated




################################################################################
################################################################################
## Node: Run

def test_run__order():
    cfg_mock  = flexmock(temp_root = "/tmp")
    node_mock = flexmock(Node())
    node_mock.should_receive("_setup").with_args(cfg_mock, "/tmp/xTMPx").ordered.once
    node_mock.should_receive("_run").with_args(cfg_mock, "/tmp/xTMPx").ordered.once
    node_mock.should_receive("_teardown").with_args(cfg_mock, "/tmp/xTMPx").ordered.once

    with monkeypatch_tmp_dir():
        node_mock.run(cfg_mock) # pylint: disable=E1103

# Ensure that the same temp dir is passed to all _ functions
def test_run__temp_dirs():
    def assert_dir(_, path):
        assert_equal(path, "/tmp/xTMPx")

    cfg_mock  = flexmock(temp_root = "/tmp")
    node_mock = flexmock(Node(),
                         _setup    = assert_dir,
                         _run      = assert_dir,
                         _teardown = assert_dir)

    with monkeypatch_tmp_dir():
        node_mock.run(cfg_mock) # pylint: disable=E1103


def test_run__exceptions():
    cfg_mock = flexmock(temp_root = "/tmp")
    def build_tests(key, exception, expectation):
        @nose.tools.raises(expectation)
        def test_function():
            node_mock = flexmock(Node())
            node_mock.should_receive(key).and_raise(exception).once
            with monkeypatch_tmp_dir():
                node_mock.run(cfg_mock) # pylint: disable=E1103

        return test_function

    for key in ('_setup', '_run', '_teardown'):
        yield build_tests(key, TypeError("The castle AAARGH!"), NodeUnhandledException)
        yield build_tests(key, NodeError("He's a very naughty boy!"), NodeError)


def test_run__error_log__node_error():
    @with_temp_folder
    def _do_test_run__error_log__node_error(temp_folder, exception):
        cfg_mock = flexmock(temp_root = temp_folder)
        node_mock = flexmock(Node())
        node_mock.should_receive("_run").and_raise(exception).once

        try:
            os.mkdir(os.path.join(temp_folder, "xTMPx"))
            with monkeypatch_tmp_dir(root = temp_folder, subfolder = "xTMPx"):
                # pylint: disable=E1103
                node_mock.run(cfg_mock) # pragma: no coverage
        except NodeError:
            log_file = os.path.join(temp_folder, "xTMPx", "pipe.errors")
            assert os.path.exists(log_file)
            assert_in("Errors =", get_file_contents(log_file))
            return
        assert False # pragma: no coverage
    yield _do_test_run__error_log__node_error, NodeError("ARGH!")
    yield _do_test_run__error_log__node_error, OSError("ARGH!")




################################################################################
################################################################################
## Node: _setup / _teardown

def test__setup__input_files():
    def _do_test__setup__input_files_exist(kwargs):
        Node(**kwargs)._setup(None, None)
    yield _do_test__setup__input_files_exist, {"executables"     : ("ls", "sh")}
    yield _do_test__setup__input_files_exist, {"input_files"     : _IN_FILES}
    yield _do_test__setup__input_files_exist, {"auxiliary_files" : _IN_FILES}

def test__setup__input_files_missing():
    @nose.tools.raises(NodeError)
    def _do_test__setup__input_files_exist(kwargs):
        Node(**kwargs)._setup(None, None)
    yield _do_test__setup__input_files_exist, {"executables"     : ("ls", "shxxxxxxxxxx")}
    yield _do_test__setup__input_files_exist, {"input_files"     : _OUT_FILES}
    yield _do_test__setup__input_files_exist, {"auxiliary_files" : _OUT_FILES}

def test__teardown__output_files():
    Node(output_files = _IN_FILES)._teardown(None, None)

@nose.tools.raises(NodeError)
def test__teardown__output_files_missing():
    Node(output_files = _OUT_FILES)._teardown(None, None)




################################################################################
################################################################################
## CommandNode: Constructor

_SIMPLE_DEPS = Node()
_SIMPLE_SUBS = Node()
_SIMPLE_CMD_MOCK = flexmock(input_files  = _IN_FILES,
                            output_files = _OUT_FILES,
                            executables  = _EXEC_FILES,
                            auxiliary_files = _AUX_FILES,
                            requirements = _REQUIREMENTS)
_SIMPLE_CMD_NODE = CommandNode(command      = _SIMPLE_CMD_MOCK,
                               subnodes     = _SIMPLE_SUBS,
                               dependencies = _SIMPLE_DEPS)

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

def test_commandnode_constructor__subnodes():
    assert_equal(_SIMPLE_CMD_NODE.subnodes, frozenset([_SIMPLE_SUBS]))

def test_commandnode_constructor__subnodes__default():
    cmd_mock = CommandNode(command = _SIMPLE_CMD_MOCK, dependencies = _SIMPLE_DEPS)
    assert_equal(cmd_mock.subnodes, frozenset())

def test_commandnode_constructor__dependencies():
    assert_equal(_SIMPLE_CMD_NODE.dependencies, frozenset([_SIMPLE_DEPS]))

def test_commandnode_constructor__dependencies__default():
    cmd_mock = CommandNode(command = _SIMPLE_CMD_MOCK, subnodes = _SIMPLE_SUBS)
    assert_equal(cmd_mock.dependencies, frozenset())



################################################################################
################################################################################
## CommandNode: run

def test_command_node__run():
    cfg_mock  = flexmock(temp_root = "/tmp")
    cmd_mock  = _build_cmd_mock()
    node_mock = flexmock(CommandNode(cmd_mock))
    node_mock.should_receive("_setup").with_args(cfg_mock, str).ordered.once
    node_mock.should_receive("_run").with_args(cfg_mock, str).ordered.once
    node_mock.should_receive("_teardown").with_args(cfg_mock, str).ordered.once
    with monkeypatch_tmp_dir():
        node_mock.run(cfg_mock) # pylint: disable=E1103




################################################################################
################################################################################
## CommandNode: _setup()

def test_commandnode_setup__files_exist():
    def _do_test_commandnode_setup(kwargs):
        cmd_mock = _build_cmd_mock(**kwargs)
        node = CommandNode(cmd_mock)
        node._setup(None, None)
    yield _do_test_commandnode_setup, {"executables"     : ("ls", "sh")}
    yield _do_test_commandnode_setup, {"input_files"     : _IN_FILES}
    yield _do_test_commandnode_setup, {"auxiliary_files" : _IN_FILES}


def test_commandnode_setup__files_missing():
    @nose.tools.raises(NodeError)
    def _do_test_commandnode_setup(kwargs):
        cmd_mock = _build_cmd_mock(**kwargs)
        node = CommandNode(cmd_mock)
        node._setup(None, None)
    yield _do_test_commandnode_setup, {"executables"     : ("ls", "shxxxxxxxxxxxxx")}
    yield _do_test_commandnode_setup, {"input_files"     : _OUT_FILES}
    yield _do_test_commandnode_setup, {"auxiliary_files" : _OUT_FILES}




################################################################################
################################################################################
## CommandNode: _run()

def test_commandnode_run__call_order():
    cmd_mock = _build_cmd_mock()
    cmd_mock.should_receive("run").with_args("xTMPx").ordered.once
    cmd_mock.should_receive("join").with_args().and_return((0,)).ordered.once
    node = CommandNode(cmd_mock)
    node._run(None, "xTMPx")

@nose.tools.raises(CmdNodeError)
def test_commandnode_run__exception_on_error():
    cmd_mock = _build_cmd_mock()
    cmd_mock.should_receive("run").ordered.once
    cmd_mock.should_receive("join").and_return((1,)).ordered.once
    node = CommandNode(cmd_mock)
    node._run(None, None)



################################################################################
################################################################################
## CommandNode: _teardown

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
                    OUT_STDOUT = os.path.join(destination, "foo.txt"))
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
@nose.tools.raises(CmdNodeError)
def test_commandnode_teardown__missing_files_in_temp(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    cmd = AtomicCmd(("echo", "-n", "1 2 3"),
                    OUT_BAR = os.path.join(destination, "bar.txt"),
                    OUT_STDOUT = os.path.join(destination, "foo.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    temp_files_before = set(os.listdir(temp_folder))
    dest_files_before = set(os.listdir(destination))

    try:
        node._teardown(None, temp_folder)
    except CmdNodeError:
        assert_equal(temp_files_before, set(os.listdir(temp_folder)))
        assert_equal(dest_files_before, set(os.listdir(destination)))
        raise

# Not all specified TEMP_ files exist at _teardown (allowed)
@with_temp_folder
def test_commandnode_teardown__missing_optional_files(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    cmd = AtomicCmd(("echo", "-n", "1 2 3"),
                    TEMP_OUT_BAR = "bar.txt",
                    OUT_STDOUT = os.path.join(destination, "foo.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    node._teardown(None, temp_folder)
    assert_equal(os.listdir(temp_folder), [])
    assert_equal(os.listdir(destination), ["foo.txt"])


# Not all required files were in place after commit
@with_temp_folder
@nose.tools.raises(NodeError)
def _test_commandnode_teardown__missing_files_in_dest(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)
    class _CmdMock(AtomicCmd):
        def commit(self, temp):
            AtomicCmd.commit(self, temp)
            os.remove(os.path.join(destination, "foo.txt"))

    cmd = _CmdMock(("touch", "%(OUT_FOO)s", "%(OUT_BAR)s"),
                   OUT_FOO = os.path.join(destination, "foo.txt"),
                   OUT_BAR = os.path.join(destination, "bar.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    node._teardown(None, temp_folder)


# Unexpected files were found in the temporary directory
@with_temp_folder
@nose.tools.raises(CmdNodeError)
def test_commandnode_teardown__extra_files_in_temp(temp_folder):
    destination, temp_folder = _setup_temp_folders(temp_folder)

    cmd = AtomicCmd(("echo", "-n", "1 2 3"),
                    OUT_STDOUT = os.path.join(destination, "foo.txt"))
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    node = CommandNode(cmd)
    set_file_contents(os.path.join(temp_folder, "bar.txt"), "1 2 3")
    temp_files_before = set(os.listdir(temp_folder))
    dest_files_before = set(os.listdir(destination))

    try:
        node._teardown(None, temp_folder)
    except CmdNodeError:
        assert_equal(temp_files_before, set(os.listdir(temp_folder)))
        assert_equal(dest_files_before, set(os.listdir(destination)))
        raise




################################################################################
################################################################################
## MetaNode: Dependencies / Subnodes

def test_metanode__nodes():
    subnodes = [Node(), Node()]
    dependencies = [Node(), Node()]
    node = MetaNode(subnodes = iter(subnodes),
                    dependencies = iter(dependencies))
    assert_equal(node.subnodes, frozenset(subnodes))
    assert_equal(node.dependencies, frozenset(dependencies))




################################################################################
################################################################################
## MetaNode: Not implemented

def test_metanode__properties():
    node = MetaNode()
    assert_equal(node.input_files, frozenset())
    assert_equal(node.output_files, frozenset())
    assert_equal(node.executables, frozenset())
    assert_equal(node.auxiliary_files, frozenset())
    assert_equal(node.requirements, frozenset())

@nose.tools.raises(MetaNodeError)
def test_metanode__is_done():
    MetaNode().is_done

@nose.tools.raises(MetaNodeError)
def test_metanode__is_outdated():
    MetaNode().is_outdated

@nose.tools.raises(MetaNodeError)
def test_metanode__run():
    MetaNode().run(None)
