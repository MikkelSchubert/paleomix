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
# pylint: disable=W0212
# Disable warnings caused by flexmock setups ("X is assigned to nothing")
# pylint: disable=W0106
import os

import nose.tools
from nose.tools import assert_equal # pylint: disable=E0611
from flexmock import flexmock
from tests.common.utils import monkeypatch, with_temp_folder

from pypeline.node import Node, MetaNode, CommandNode, NodeError, NodeUnhandledException



class monkeypatch_tmp_dir:
    def __init__(self, root = "/tmp", subfolder = "xTMPx", require_rmdir = True):
        self._mkdir_patch = monkeypatch("pypeline.common.fileutils.create_temp_dir", self._mkdir)
        self._rmdir_patch = monkeypatch("os.rmdir", self._rmdir)
        self._root_dir = root
        self._sub_dir  = subfolder
        self._mkdir_called = False
        self._rmdir_called = False
        self._require_rmdir = require_rmdir

    def __enter__(self):
        self._mkdir_patch.__enter__()
        self._rmdir_patch.__enter__()
        return self

    def __exit__(self, type, value, traceback):
        self._mkdir_patch.__exit__(type, value, traceback)
        self._rmdir_patch.__exit__(type, value, traceback)
        assert self._mkdir_called
        assert self._rmdir_called or not self._require_rmdir

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
_AUX_FILES   = frozenset(("tests/data/rCRS.fasta.fai",))
_REQUIREMENTS = frozenset((id,))


def _build_cmd_mock(input_files  = (), output_files = (), executables  = (), auxiliary_files = (), requirements = (),
                    expected_temp_files = (), optional_temp_files = ()):
    return flexmock(input_files     = frozenset(input_files),
                    output_files    = frozenset(output_files),
                    executables     = frozenset(executables),
                    auxiliary_files = frozenset(auxiliary_files),
                    requirements    = frozenset(requirements),
                    expected_temp_files = frozenset(expected_temp_files),
                    optional_temp_files = frozenset(optional_temp_files))

################################################################################
################################################################################
## Node: Constructor tests

def test_constructor__single_input_file():
    in_file = iter(_IN_FILES).next()
    my_node = Node(input_files  = in_file)
    assert_equal(my_node.input_files, frozenset([in_file]))

def test_constructor__input_files():
    my_node = Node(input_files  = _IN_FILES)
    assert_equal(my_node.input_files, _IN_FILES)



def test_constructor__single_output_file():
    out_file = iter(_OUT_FILES).next()
    my_node = Node(output_files  = out_file)
    assert_equal(my_node.output_files, frozenset([out_file]))

def test_constructor__output_files():
    my_node = Node(output_files  = _OUT_FILES)
    assert_equal(my_node.output_files, _OUT_FILES)



def test_constructor__single_executable():
    out_file = iter(_OUT_FILES).next()
    my_node = Node(executables  = out_file)
    assert_equal(my_node.executables, frozenset([out_file]))

def test_constructor__executables():
    my_node = Node(executables  = _OUT_FILES)
    assert_equal(my_node.executables, _OUT_FILES)



def test_constructor__description():
    my_node = Node(description  = _DESCRIPTION)
    assert_equal(str(my_node), _DESCRIPTION)



def test_constructor__single_subnode():
    subnode = Node()
    my_node  = Node(subnodes = subnode)
    assert_equal(my_node.subnodes, frozenset([subnode]))

def test_constructor__mult_subnodes():
    subnodes = set([Node(), Node()])
    my_node  = Node(subnodes = subnodes)
    assert_equal(my_node.subnodes, subnodes)

@nose.tools.raises(TypeError)
def test_constructor__subnode_not_a_node():
    Node(subnodes = (1,))




def test_constructor__single_dependencies():
    dependencies = Node()
    my_node  = Node(dependencies = dependencies)
    assert_equal(my_node.dependencies, frozenset([dependencies]))

def test_constructor__mult_dependencies():
    dependencies = set([Node(), Node()])
    my_node  = Node(dependencies = dependencies)
    assert_equal(my_node.dependencies, dependencies)

@nose.tools.raises(TypeError)
def test_constructor__dependency_not_a_node():
    Node(dependencies = (1,))




################################################################################
################################################################################
## Node: Properties tests

def test_is_done__no_output():
    assert Node().is_done

def test_is_done__output_missing():
    my_node = Node(output_files = "tests/data/missing_file")
    assert not my_node.is_done

def test_is_done__output_exists():
    my_node = Node(output_files = "tests/data/empty_file_1")
    assert my_node.is_done

def test_is_done__updates():
    my_node = Node(output_files = "tests/data/missing_file")
    assert not my_node.is_done
    my_node.output_files = ("tests/data/empty_file_1",)
    assert my_node.is_done

def test_is_done__subnodes_not_done():
    my_subnode = Node(output_files = "tests/data/missing_file")
    my_node    = Node(subnodes = my_subnode)
    assert not my_node.is_done

def test_is_done__subnodes_is_done():
    my_subnode = Node(output_files = "tests/data/empty_file_1")
    my_node    = Node(subnodes = my_subnode)
    assert my_node.is_done

def test_is_done__subnodes_updated():
    my_subnode = Node(output_files = "tests/data/missing_file")
    my_node    = Node(subnodes = my_subnode)
    assert not my_node.is_done
    my_subnode.output_files = ("tests/data/empty_file_1",)
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




def test_run__order():
    cfg_mock  = flexmock(temp_root = "/tmp")
    node_mock = flexmock(Node())
    node_mock.should_receive("_setup").with_args(cfg_mock, "/tmp/xTMPx").ordered.once
    node_mock.should_receive("_run").with_args(cfg_mock, "/tmp/xTMPx").ordered.once
    node_mock.should_receive("_teardown").with_args(cfg_mock, "/tmp/xTMPx").ordered.once

    with monkeypatch_tmp_dir():
        node_mock.run(cfg_mock) # pylint: disable=E1103

def test_run__temp_dirs():
    paths = []
    def assert_dir(_, path):
        paths.append(path)

    cfg_mock  = flexmock(temp_root = "/tmp")
    node_mock = flexmock(Node(),
                         _setup    = assert_dir,
                         _run      = assert_dir,
                         _teardown = assert_dir)

    with monkeypatch_tmp_dir():
        node_mock.run(cfg_mock) # pylint: disable=E1103
    assert len(set(paths)) == 1, paths


def test_run__exceptions():
    cfg_mock = flexmock(temp_root = "/tmp")
    def build_tests(key, exception, expectation):
        @nose.tools.raises(expectation)
        def test_function():
            node_mock = flexmock(Node())
            node_mock.should_receive(key).and_raise(exception).once
            with monkeypatch_tmp_dir(require_rmdir = False):
                node_mock.run(cfg_mock) # pylint: disable=E1103

        return test_function

    for key in ('_setup', '_run', '_teardown'):
        yield build_tests(key, TypeError("The castle AAARGH!"), NodeUnhandledException)
        yield build_tests(key, NodeError("He's a very naughty boy!"), NodeError)


def test__setup__input_files_exist():
    Node(input_files  = _IN_FILES)._setup(None, None)

@nose.tools.raises(NodeError)
def test__setup__input_files_missing():
    Node(input_files  = _OUT_FILES)._setup(None, None)


def test__setup__output_files_exist():
    Node(output_files  = _IN_FILES)._teardown(None, None)

@nose.tools.raises(NodeError)
def test__setup__output_files_missing():
    Node(output_files  = _OUT_FILES)._teardown(None, None)




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
                               description  = "SimpleCommand",
                               subnodes     = _SIMPLE_SUBS,
                               dependencies = _SIMPLE_DEPS)

def test_commandnode_constructor__description():
    assert_equal(str(_SIMPLE_CMD_NODE), "SimpleCommand")

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

def test_commandnode_constructor__dependencies():
    assert_equal(_SIMPLE_CMD_NODE.dependencies, frozenset([_SIMPLE_DEPS]))


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

def test_commandnode_setup__executable_missing():
    cmd_mock = _build_cmd_mock(executables  = ("ls", "sh"))
    node = CommandNode(cmd_mock)
    node._setup(None, None)

@nose.tools.raises(NodeError)
def test_commandnode_setup__executables_exist():
    cmd_mock = _build_cmd_mock(executables  = ("ls", "shxxxxxxxxxxxxx"))
    node = CommandNode(cmd_mock)
    node._setup(None, None)

def test_commandnode_setup__input_exist():
    cmd_mock = _build_cmd_mock(input_files  = _IN_FILES)
    node = CommandNode(cmd_mock)
    node._setup(None, None)

@nose.tools.raises(NodeError)
def test_commandnode_setup__input_missing():
    cmd_mock = _build_cmd_mock(input_files  = _OUT_FILES)
    node = CommandNode(cmd_mock)
    node._setup(None, None)




################################################################################
################################################################################
## CommandNode: _run()

def test_commandnode_run__call_order():
    cmd_mock = _build_cmd_mock()
    cmd_mock.should_receive("run").with_args("xTMPx").ordered.once
    cmd_mock.should_receive("join").with_args().and_return((0,)).ordered.once
    node = CommandNode(cmd_mock)
    node._run(None, "xTMPx")

@nose.tools.raises(NodeError)
def test_commandnode_run__exception_on_error():
    cmd_mock = _build_cmd_mock()
    cmd_mock.should_receive("run").ordered.once
    cmd_mock.should_receive("join").and_return((1,)).ordered.once
    node = CommandNode(cmd_mock)
    node._run(None, None)



################################################################################
################################################################################
## CommandNode: _teardown

@with_temp_folder
def test_commandnode_teardown__output_exists(temp_folder):
    cmd_mock = _build_cmd_mock(output_files = _IN_FILES)
    cmd_mock.should_receive("commit").once
    node = CommandNode(cmd_mock)
    node._teardown(None, temp_folder)

@with_temp_folder
@nose.tools.raises(NodeError)
def test_commandnode_teardown__output_missing(temp_folder):
    cmd_mock = _build_cmd_mock(output_files = _OUT_FILES)
    cmd_mock.should_receive("commit").once
    node = CommandNode(cmd_mock)
    node._teardown(None, temp_folder)
