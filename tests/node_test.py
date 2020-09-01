#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import logging
import os
import random
import uuid

from unittest.mock import call, Mock

import pytest

import paleomix.node

from paleomix.atomiccmd.command import AtomicCmd
from paleomix.node import (
    Node,
    CommandNode,
    NodeError,
    NodeUnhandledException,
    CmdNodeError,
)
from paleomix.common.utilities import safe_coerce_to_frozenset


def test_dir():
    return os.path.dirname(__file__)


def test_file(*args):
    return os.path.join(test_dir(), "data", *args)


def first(values):
    return random.choice(tuple(values))


def _CommandNodeWrap(**kwargs):
    return CommandNode(command=AtomicCmd("true"), **kwargs)


_NODE_TYPES = (Node, _CommandNodeWrap)


_DESCRIPTION = "My description of a node"
_IN_FILES = frozenset((test_file("empty_file_1"), test_file("empty_file_2")))
_OUT_FILES = frozenset(
    (test_file("missing_out_file_1"), test_file("missing_out_file_2"))
)
_EXEC_FILES = frozenset(("ls", "sh"))
_AUX_FILES = frozenset((test_file("rCRS.fasta"), test_file("rCRS.fasta.fai")))
_REQUIREMENTS = frozenset((id, str))
_EMPTY_FILE = test_file("empty_file_1")


def _build_cmd_mock(
    input_files=_IN_FILES,
    output_files=(),
    executables=(),
    auxiliary_files=(),
    requirements=(),
    optional_temp_files=(),
    return_codes=(0,),
):
    cmd = Mock(
        input_files=frozenset(input_files),
        output_files=frozenset(output_files),
        executables=frozenset(executables),
        auxiliary_files=frozenset(auxiliary_files),
        requirements=frozenset(requirements),
        expected_temp_files=frozenset(map(os.path.basename, output_files)),
        optional_temp_files=frozenset(optional_temp_files),
    )
    cmd.join.return_value = return_codes

    return cmd


###############################################################################
###############################################################################
# Node: Constructor: File sets


_CONSTUCTOR_SINGLE_VALUES = (
    # Single values
    ("input_files", first(_IN_FILES)),
    ("output_files", first(_OUT_FILES)),
    ("executables", first(_EXEC_FILES)),
    ("auxiliary_files", first(_AUX_FILES)),
    # Single value in list
    ("input_files", [first(_IN_FILES)]),
    ("output_files", [first(_OUT_FILES)]),
    ("executables", [first(_EXEC_FILES)]),
    ("auxiliary_files", [first(_AUX_FILES)]),
    # Multiple values in list
    ("input_files", _IN_FILES),
    ("output_files", _OUT_FILES),
    ("executables", _EXEC_FILES),
    ("auxiliary_files", _AUX_FILES),
)


@pytest.mark.parametrize("key, value", _CONSTUCTOR_SINGLE_VALUES)
def test_constructor(key, value):
    defaults = {"input_files": _EMPTY_FILE}
    defaults[key] = value
    node = Node(**defaults)
    expected = safe_coerce_to_frozenset(value)
    assert getattr(node, key) == expected


_CONSTUCTOR_INVALID_VALUES = (
    ("input_files", [id]),
    ("output_files", [-1]),
    ("executables", [{}]),
    ("auxiliary_files", [1.3]),
)


@pytest.mark.parametrize("key, value", _CONSTUCTOR_INVALID_VALUES)
def test_constructor__invalid_values(key, value):
    with pytest.raises(TypeError):
        Node(**{key: value})


###############################################################################
###############################################################################
# Node: Constructor: Requirements


def test_constructor__requirements():
    node = Node(requirements=id)
    assert node.requirements == frozenset([id])
    node = Node(requirements=[id])
    assert node.requirements == frozenset([id])
    node = Node(requirements=[id, str])
    assert node.requirements == frozenset([id, str])


@pytest.mark.parametrize("value", (17, {}, "867-5309"))
def test_constructor__requirements__wrong_type(value):
    with pytest.raises(TypeError):
        Node(requirements=value)


###############################################################################
###############################################################################
# Node: Constructor: Dependencies


def test_constructor__nodes_is_none():
    my_node = Node(dependencies=None)
    assert my_node.dependencies == frozenset()


def test_constructor__single_node():
    sub_node = Node()
    my_node = Node(dependencies=sub_node)
    assert my_node.dependencies == frozenset([sub_node])


def test_constructor__iterable():
    sub_nodes = [Node(), Node()]
    my_node = Node(dependencies=iter(sub_nodes))
    assert my_node.dependencies == frozenset(sub_nodes)


def test_constructor__not_a_node():
    with pytest.raises(TypeError):
        Node(dependencies=(1,))


###############################################################################
###############################################################################
# *Node: Description


@pytest.mark.parametrize("cls", _NODE_TYPES)
def test_constructor__description(cls):
    my_node = cls(description=_DESCRIPTION)
    assert str(my_node) == _DESCRIPTION


@pytest.mark.parametrize("cls", _NODE_TYPES)
def test_constructor__description__default(cls):
    my_node = cls()
    assert str(my_node) == repr(my_node)


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("value", (1, {}))
def test_constructor__description__non_string(cls, value):
    with pytest.raises(TypeError):
        cls(description=value)


###############################################################################
###############################################################################
# *Node: Constructor tests: #threads


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("nthreads", (1, 3))
def test_constructor__threads(cls, nthreads):
    node = cls(threads=nthreads)
    assert node.threads == nthreads


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("nthreads", (-1, 0))
def test_constructor__threads_invalid_range(cls, nthreads):
    with pytest.raises(ValueError):
        cls(threads=nthreads)


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("nthreads", ("1", {}, 2.7))
def test_constructor__threads_invalid_type(cls, nthreads):
    with pytest.raises(TypeError):
        cls(threads=nthreads)


###############################################################################
###############################################################################
# Node: Run

_DUMMY_TEMP_ROOT = "/xyz/tmp"
_DUMMY_TEMP = os.path.join(_DUMMY_TEMP_ROOT, "xTMPx")


def test_run__order():
    cfg_mock = Mock(temp_root=_DUMMY_TEMP_ROOT)
    node_mock = Mock()

    node = Node()
    node._create_temp_dir = node_mock._create_temp_dir
    node._create_temp_dir.return_value = _DUMMY_TEMP
    node._setup = node_mock._setup
    node._run = node_mock._run
    node._teardown = node_mock._teardown
    node._remove_temp_dir = node_mock._remove_temp_dir

    node.run(cfg_mock)

    node_mock.mock_calls == [
        call._create_temp_dir(cfg_mock),
        call._setup(cfg_mock, _DUMMY_TEMP),
        call._run(cfg_mock, _DUMMY_TEMP),
        call._teardown(cfg_mock, _DUMMY_TEMP),
        call._remove_temp_dir(_DUMMY_TEMP),
    ]


_EXCEPTIONS = (
    (TypeError("The castle AAARGH!"), NodeUnhandledException),
    (NodeError("He's a very naughty boy!"), NodeError),
)


@pytest.mark.parametrize("key", ("_setup", "_run", "_teardown"))
@pytest.mark.parametrize("exception, expectation", _EXCEPTIONS)
def test_run__exceptions(key, exception, expectation):
    mock = Mock()
    node = Node()
    node._create_temp_dir = mock._create_temp_dir
    node._create_temp_dir.return_value = _DUMMY_TEMP

    setattr(node, key, getattr(mock, key))
    getattr(node, key).side_effect = exception

    cfg_mock = Mock(temp_root=_DUMMY_TEMP_ROOT)
    with pytest.raises(expectation):
        node.run(cfg_mock)

    assert mock.mock_calls == [
        call._create_temp_dir(cfg_mock),
        getattr(call, key)(cfg_mock, _DUMMY_TEMP),
    ]


def test_run__exception__create_temp_dir():
    cfg_mock = Mock(temp_root=_DUMMY_TEMP_ROOT)
    node_mock = Node()
    node_mock._create_temp_dir = Mock()
    node_mock._create_temp_dir.side_effect = OSError()

    with pytest.raises(NodeUnhandledException):
        node_mock.run(cfg_mock)
    assert node_mock._create_temp_dir.mock_calls == [call(cfg_mock)]


def test_run__exception__remove_temp_dir():
    cfg_mock = Mock(temp_root=_DUMMY_TEMP_ROOT)
    mock = Mock()
    node_mock = Node()
    node_mock._create_temp_dir = mock._create_temp_dir
    node_mock._create_temp_dir.return_value = _DUMMY_TEMP
    node_mock._remove_temp_dir = mock._remove_temp_dir
    node_mock._remove_temp_dir.side_effect = OSError()

    with pytest.raises(NodeUnhandledException):
        node_mock.run(cfg_mock)
    assert mock.mock_calls == [
        call._create_temp_dir(cfg_mock),
        call._remove_temp_dir(_DUMMY_TEMP),
    ]


@pytest.mark.parametrize("exception", (NodeError, OSError))
def test_run__error_log__node_error(tmp_path, exception):
    temp = tmp_path / "xTMPx"
    mock = Mock()
    cfg_mock = Mock(temp_root=tmp_path)
    node_mock = Node()
    node_mock._create_temp_dir = mock._create_temp_dir
    node_mock._create_temp_dir.return_value = str(temp)
    node_mock._run = mock._run
    node_mock._run.side_effect = exception("ARGH!")

    temp.mkdir()
    with pytest.raises(NodeError):
        node_mock.run(cfg_mock)
    log_file = tmp_path / "xTMPx" / "pipe.errors"
    assert log_file.exists()

    error_text = log_file.read_text()
    assert "Errors =" in error_text
    assert paleomix.__version__ in error_text

    assert mock.mock_calls == [
        call._create_temp_dir(cfg_mock),
        call._run(cfg_mock, str(temp)),
    ]


###############################################################################
###############################################################################
# Node: _setup / _teardown

_INPUT_FILES_EXIST = (
    {"executables": ("ls", "sh")},
    {"input_files": _IN_FILES},
    {"auxiliary_files": _IN_FILES},
)


@pytest.mark.parametrize("kwargs", _INPUT_FILES_EXIST)
def test__setup__input_files(kwargs):
    Node(**kwargs)._setup(None, None)


_INPUT_FILES_MISSING = (
    {"executables": ("ls", "shxxxx")},
    {"input_files": _OUT_FILES},
    {"auxiliary_files": _OUT_FILES},
)


@pytest.mark.parametrize("kwargs", _INPUT_FILES_MISSING)
def test__setup__input_files_missing(kwargs):
    with pytest.raises(NodeError):
        Node(**kwargs)._setup(None, None)


def test__teardown__output_files():
    Node(input_files=_EMPTY_FILE, output_files=_IN_FILES)._teardown(None, None)


def test__teardown__output_files_missing():
    node = Node(input_files=_EMPTY_FILE, output_files=_OUT_FILES)
    with pytest.raises(NodeError):
        node._teardown(None, None)


###############################################################################
# Node._remove_temp_dir


def test_node_remove_temp_dir__empty_dir(tmp_path, caplog):
    with caplog.at_level(logging.WARNING):
        node = Node()
        node._remove_temp_dir(tmp_path)

        assert not tmp_path.exists()
        assert not caplog.messages


def test_node_remove_temp_dir__extranous_files(tmp_path, caplog):
    tmp_file = tmp_path / str(uuid.uuid4())
    tmp_file.touch()

    with caplog.at_level(logging.WARNING):
        node = Node()
        node._remove_temp_dir(tmp_path)

        assert not tmp_path.exists()
        assert (
            paleomix.node.__name__,
            logging.WARNING,
            "Unexpected file in temporary directory: %r" % (str(tmp_file),),
        ) in caplog.record_tuples


###############################################################################
# Node._collect_files


def test_node_collect_files__empty_folder(tmp_path):
    assert list(Node._collect_files(tmp_path)) == []


def test_node_collect_files__root_files(tmp_path):
    (tmp_path / "foo.txt").touch()
    (tmp_path / "bar.txt").touch()

    assert sorted(Node._collect_files(tmp_path)) == [
        "bar.txt",
        "foo.txt",
    ]


def test_node_collect_files__files_and_folders(tmp_path):
    (tmp_path / "foo1.txt").touch()
    (tmp_path / "bar1").mkdir()
    (tmp_path / "bar1" / "foo2.txt").touch()
    (tmp_path / "bar1" / "bar2").mkdir()
    (tmp_path / "bar1" / "bar2" / "foo3.txt").touch()
    (tmp_path / "bar1" / "bar2" / "foo4.txt").touch()

    assert sorted(Node._collect_files(tmp_path)) == [
        "bar1/bar2/foo3.txt",
        "bar1/bar2/foo4.txt",
        "bar1/foo2.txt",
        "foo1.txt",
    ]


###############################################################################
###############################################################################
# CommandNode: Constructor

_SIMPLE_DEPS = Node()
_SIMPLE_SUBS = Node()
_SIMPLE_CMD_MOCK = Mock(
    input_files=_IN_FILES,
    output_files=_OUT_FILES,
    executables=_EXEC_FILES,
    auxiliary_files=_AUX_FILES,
    requirements=_REQUIREMENTS,
)
_SIMPLE_CMD_NODE = CommandNode(command=_SIMPLE_CMD_MOCK, dependencies=_SIMPLE_DEPS)


def test_commandnode_constructor__input_files():
    assert _SIMPLE_CMD_NODE.input_files == _IN_FILES


def test_commandnode_constructor__output_files():
    assert _SIMPLE_CMD_NODE.output_files == _OUT_FILES


def test_commandnode_constructor__auxiliary_files():
    assert _SIMPLE_CMD_NODE.auxiliary_files == _AUX_FILES


def test_commandnode_constructor__executables():
    assert _SIMPLE_CMD_NODE.executables == _EXEC_FILES


def test_commandnode_constructor__requirements():
    assert _SIMPLE_CMD_NODE.requirements == _REQUIREMENTS


def test_commandnode_constructor__dependencies():
    assert _SIMPLE_CMD_NODE.dependencies == frozenset([_SIMPLE_DEPS])


def test_commandnode_constructor__dependencies__default():
    cmd_mock = CommandNode(command=_SIMPLE_CMD_MOCK)
    assert cmd_mock.dependencies == frozenset()


###############################################################################
###############################################################################
# CommandNode: run


def test_command_node__run():
    cfg_mock = Mock(temp_root=_DUMMY_TEMP_ROOT)
    mock = _build_cmd_mock()

    node_mock = CommandNode(mock)
    node_mock._create_temp_dir = mock._test_node_._create_temp_dir
    node_mock._create_temp_dir.return_value = _DUMMY_TEMP
    node_mock._setup = mock._test_node_._setup
    node_mock._teardown = mock._test_node_._teardown
    node_mock._remove_temp_dir = mock._test_node_._remove_temp_dir

    node_mock.run(cfg_mock)

    assert mock.mock_calls == [
        call._test_node_._create_temp_dir(cfg_mock),
        call._test_node_._setup(cfg_mock, _DUMMY_TEMP),
        call.run(_DUMMY_TEMP),
        call.join(),
        call._test_node_._teardown(cfg_mock, _DUMMY_TEMP),
        call._test_node_._remove_temp_dir(_DUMMY_TEMP),
    ]


###############################################################################
###############################################################################
# CommandNode: _setup()

_SETUP_FILES_EXIST = (
    {"executables": ("ls", "sh")},
    {"input_files": _IN_FILES},
    {"auxiliary_files": _IN_FILES},
)


@pytest.mark.parametrize("kwargs", _INPUT_FILES_EXIST)
def test_commandnode_setup__files_exist(kwargs):
    cmd_mock = _build_cmd_mock(**kwargs)
    node = CommandNode(cmd_mock)
    node._setup(None, None)


_SETUP_FILES_MISSING = (
    {"executables": ("ls", "shxxxxxxxxxxx")},
    {"input_files": _OUT_FILES},
    {"auxiliary_files": _OUT_FILES},
)


@pytest.mark.parametrize("kwargs", _INPUT_FILES_MISSING)
def test_commandnode_setup__files_missing(kwargs):
    cmd_mock = _build_cmd_mock(**kwargs)
    node = CommandNode(cmd_mock)
    with pytest.raises(NodeError):
        node._setup(None, None)


###############################################################################
###############################################################################
# CommandNode: _run()


def test_commandnode_run__call_order():
    cmd_mock = _build_cmd_mock()
    node = CommandNode(cmd_mock)
    node._run(None, "xTMPx")

    assert cmd_mock.mock_calls == [call.run("xTMPx"), call.join()]


def test_commandnode_run__exception_on_error():
    cmd_mock = _build_cmd_mock(return_codes=(1,))
    node = CommandNode(cmd_mock)
    with pytest.raises(CmdNodeError):
        node._run(None, "xTMPx")

    assert cmd_mock.mock_calls == [call.run("xTMPx"), call.join()]


###############################################################################
###############################################################################
# CommandNode: _teardown


def _setup_temp_folders(tmp_path):
    destination = tmp_path / "dst"
    tmp_path = tmp_path / "tmp"
    tmp_path.mkdir(parents=True, exist_ok=True)
    destination.mkdir(parents=True, exist_ok=True)
    return destination, tmp_path


# Commit is called on the command obj
def test_commandnode_teardown__commit(tmp_path):
    cmd_mock = _build_cmd_mock()
    node = CommandNode(cmd_mock)
    node._teardown(None, tmp_path)
    assert cmd_mock.mock_calls == [call.commit(tmp_path)]


# Files exist in temp folder, and in destination after commit
def test_commandnode_teardown(tmp_path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    cmd = AtomicCmd(
        ("echo", "-n", "1 2 3"),
        IN_DUMMY=_EMPTY_FILE,
        OUT_STDOUT=str(destination / "foo.txt"),
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    assert (tmp_path / "foo.txt").exists()
    assert not (destination / "foo.txt").exists()
    node._teardown(None, tmp_path)
    assert not (tmp_path / "foo.txt").exists()
    assert (destination / "foo.txt").exists()


# Not all required files have been generated (atomic)
def test_commandnode_teardown__missing_files_in_temp(tmp_path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    cmd = AtomicCmd(
        ("echo", "-n", "1 2 3"),
        IN_DUMMY=_EMPTY_FILE,
        OUT_BAR=str(destination / "bar.txt"),
        OUT_STDOUT=str(destination / "foo.txt"),
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    temp_files_before = set(tmp_path.iterdir())
    dest_files_before = set(destination.iterdir())

    with pytest.raises(CmdNodeError):
        node._teardown(None, tmp_path)
    assert temp_files_before == set(tmp_path.iterdir())
    assert dest_files_before == set(destination.iterdir())


# Not all specified TEMP_ files exist at _teardown (allowed)
def test_commandnode_teardown__missing_optional_files(tmp_path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    cmd = AtomicCmd(
        ("echo", "-n", "1 2 3"),
        IN_DUMMY=_EMPTY_FILE,
        TEMP_OUT_BAR="bar.txt",
        OUT_STDOUT=str(destination / "foo.txt"),
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    node._teardown(None, tmp_path)
    assert os.listdir(str(tmp_path)) == []
    assert os.listdir(str(destination)) == ["foo.txt"]


# Not all required files were in place after commit
def test_commandnode_teardown__missing_files_in_dest(tmp_path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    class _CmdMock(AtomicCmd):
        def commit(self, temp):
            AtomicCmd.commit(self, temp)
            (destination / "foo.txt").unlink()

    cmd = _CmdMock(
        ("touch", "%(OUT_FOO)s", "%(OUT_BAR)s"),
        IN_DUMMY=_EMPTY_FILE,
        OUT_FOO=str(destination / "foo.txt"),
        OUT_BAR=str(destination / "bar.txt"),
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    with pytest.raises(NodeError):
        node._teardown(None, tmp_path)


def test_commandnode_teardown__extra_files_in_temp(tmp_path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    unexpected_file = tmp_path / "unexpected_file.txt"
    unexpected_file.write_text("1 2 3")

    cmd = AtomicCmd(
        ("echo", "1 2 3"),
        IN_DUMMY=_EMPTY_FILE,
        OUT_STDOUT=str(destination / "foo.txt"),
    )

    node = CommandNode(cmd)
    node._run(None, tmp_path)
    node._teardown(None, tmp_path)

    assert list(tmp_path.iterdir()) == [unexpected_file]
    assert list(destination.iterdir()) == [destination / "foo.txt"]
