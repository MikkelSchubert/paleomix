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
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence, Tuple, Type, TypeVar
from unittest.mock import Mock, call

import pytest

import paleomix.node
from paleomix.common.command import AtomicCmd, InputFile, OutputFile, TempOutputFile
from paleomix.common.utilities import safe_coerce_to_frozenset
from paleomix.common.versions import Requirement
from paleomix.node import (
    CmdNodeError,
    CommandNode,
    Node,
    NodeError,
    NodeUnhandledException,
)

T = TypeVar("T")


def _test_dir():
    return os.path.dirname(__file__)


def _test_file(*args: str):
    return os.path.join(_test_dir(), "data", *args)


def choice(values: Iterable[T]) -> T:
    return random.choice(tuple(values))


class CommandNodeWithDefaultCommand(CommandNode):
    def __init__(
        self,
        description: Optional[str] = None,
        threads: int = 1,
        dependencies: Iterable[Node] = (),
    ):
        super().__init__(AtomicCmd("true"), description, threads, dependencies)


_NODE_TYPES: Tuple[Type[Node], ...] = (Node, CommandNodeWithDefaultCommand)


_DESCRIPTION = "My description of a node"
_IN_FILES = frozenset((_test_file("empty_file_1"), _test_file("empty_file_2")))
_OUT_FILES = frozenset(
    (_test_file("missing_out_file_1"), _test_file("missing_out_file_2"))
)
_EXEC_FILES = frozenset(("ls", "sh"))
_AUX_FILES = frozenset((_test_file("rCRS.fasta"), _test_file("rCRS.fasta.fai")))
_REQUIREMENT_1 = Requirement(["bwa"], r"", "")
_REQUIREMENT_2 = Requirement(["bowtie"], r"", "")
_REQUIREMENTS = frozenset((_REQUIREMENT_1, _REQUIREMENT_2))
_EMPTY_FILE = _test_file("empty_file_1")


def _build_cmd_mock(
    input_files: Iterable[str] = _IN_FILES,
    output_files: Iterable[str] = (),
    executables: Iterable[str] = (),
    auxiliary_files: Iterable[str] = (),
    requirements: Iterable[str] = (),
    optional_temp_files: Iterable[str] = (),
    return_codes: Sequence[int] = (0,),
):
    cmd = Mock(
        input_files=frozenset(input_files),
        output_files=frozenset(output_files),
        executables=frozenset(executables),
        auxiliary_files=frozenset(auxiliary_files),
        requirements=frozenset(requirements),
        expected_temp_files=frozenset(os.path.basename(f) for f in output_files),
        optional_temp_files=frozenset(optional_temp_files),
    )
    cmd.join.return_value = return_codes

    return cmd


###############################################################################
###############################################################################
# Node: Constructor: File sets


_CONSTUCTOR_SINGLE_VALUES = (
    # Single values
    ("input_files", choice(_IN_FILES)),
    ("output_files", choice(_OUT_FILES)),
    ("executables", choice(_EXEC_FILES)),
    ("auxiliary_files", choice(_AUX_FILES)),
    # Single value in list
    ("input_files", [choice(_IN_FILES)]),
    ("output_files", [choice(_OUT_FILES)]),
    ("executables", [choice(_EXEC_FILES)]),
    ("auxiliary_files", [choice(_AUX_FILES)]),
    # Multiple values in list
    ("input_files", _IN_FILES),
    ("output_files", _OUT_FILES),
    ("executables", _EXEC_FILES),
    ("auxiliary_files", _AUX_FILES),
)


@pytest.mark.parametrize("key, value", _CONSTUCTOR_SINGLE_VALUES)
def test_constructor(key: str, value: str):
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
def test_constructor__invalid_values(key: str, value: Any):
    with pytest.raises(TypeError):
        Node(**{key: value})


###############################################################################
###############################################################################
# Node: Constructor: Requirements


def test_constructor__requirements():
    node = Node(requirements=_REQUIREMENT_1)
    assert node.requirements == frozenset([_REQUIREMENT_1])
    node = Node(requirements=[_REQUIREMENT_1])
    assert node.requirements == frozenset([_REQUIREMENT_1])
    node = Node(requirements=[_REQUIREMENT_1, _REQUIREMENT_2])
    assert node.requirements == frozenset([_REQUIREMENT_1, _REQUIREMENT_2])


@pytest.mark.parametrize("value", (17, "867-5309"))
def test_constructor__requirements__wrong_type(value: Any):
    with pytest.raises(TypeError):
        Node(requirements=value)  # type: ignore


###############################################################################
###############################################################################
# Node: Constructor: Dependencies


def test_constructor__nodes_is_none():
    with pytest.raises(TypeError):
        Node(dependencies=None)  # type: ignore


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
        Node(dependencies=(1,))  # type: ignore


###############################################################################
###############################################################################
# *Node: Description


@pytest.mark.parametrize("cls", _NODE_TYPES)
def test_constructor__description(cls: Type[Node]):
    my_node = cls(description=_DESCRIPTION)
    assert str(my_node) == _DESCRIPTION


@pytest.mark.parametrize("cls", _NODE_TYPES)
def test_constructor__description__default(cls: Type[Node]):
    my_node = cls()
    assert str(my_node) == repr(my_node)


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("value", (1, {}))
def test_constructor__description__non_string(cls: Type[Node], value: Any):
    with pytest.raises(TypeError):
        cls(description=value)  # type: ignore


###############################################################################
###############################################################################
# *Node: Constructor tests: #threads


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("nthreads", (1, 3))
def test_constructor__threads(cls: Type[Node], nthreads: int):
    node = cls(threads=nthreads)
    assert node.threads == nthreads


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("nthreads", (-1, 0))
def test_constructor__threads_invalid_range(cls: Type[Node], nthreads: int):
    with pytest.raises(ValueError):
        cls(threads=nthreads)


@pytest.mark.parametrize("cls", _NODE_TYPES)
@pytest.mark.parametrize("nthreads", ("1", {}, 2.7))
def test_constructor__threads_invalid_type(cls: Type[Node], nthreads: int):
    with pytest.raises(TypeError):
        cls(threads=nthreads)


###############################################################################
###############################################################################
# Node: Run

_DUMMY_TEMP_ROOT = "/xyz/tmp"
_DUMMY_TEMP = os.path.join(_DUMMY_TEMP_ROOT, "xTMPx")


def test_run__order():
    node_mock = Mock()

    node = Node()
    node._create_temp_dir = node_mock._create_temp_dir
    node._create_temp_dir.return_value = _DUMMY_TEMP
    node._setup = node_mock._setup
    node._run = node_mock._run
    node._teardown = node_mock._teardown
    node._remove_temp_dir = node_mock._remove_temp_dir

    node.run(_DUMMY_TEMP_ROOT)

    node_mock.mock_calls == [
        call._create_temp_dir(_DUMMY_TEMP_ROOT),
        call._setup(_DUMMY_TEMP),
        call._run(_DUMMY_TEMP),
        call._teardown(_DUMMY_TEMP),
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

    with pytest.raises(expectation):
        node.run(_DUMMY_TEMP_ROOT)

    assert mock.mock_calls == [
        call._create_temp_dir(_DUMMY_TEMP_ROOT),
        getattr(call, key)(_DUMMY_TEMP),
    ]


def test_run__exception__create_temp_dir():
    node_mock = Node()
    node_mock._create_temp_dir = Mock()
    node_mock._create_temp_dir.side_effect = OSError()

    with pytest.raises(NodeUnhandledException):
        node_mock.run(_DUMMY_TEMP_ROOT)
    assert node_mock._create_temp_dir.mock_calls == [call(_DUMMY_TEMP_ROOT)]


def test_run__exception__remove_temp_dir():
    mock = Mock()
    node_mock = Node()
    node_mock._create_temp_dir = mock._create_temp_dir
    node_mock._create_temp_dir.return_value = _DUMMY_TEMP
    node_mock._remove_temp_dir = mock._remove_temp_dir
    node_mock._remove_temp_dir.side_effect = OSError()

    with pytest.raises(NodeUnhandledException):
        node_mock.run(_DUMMY_TEMP_ROOT)
    assert mock.mock_calls == [
        call._create_temp_dir(_DUMMY_TEMP_ROOT),
        call._remove_temp_dir(_DUMMY_TEMP),
    ]


@pytest.mark.parametrize("exception", (NodeError, OSError))
def test_run__error_log__node_error(tmp_path, exception):
    temp = tmp_path / "xTMPx"
    mock = Mock()
    node_mock = Node()
    node_mock._create_temp_dir = mock._create_temp_dir
    node_mock._create_temp_dir.return_value = str(temp)
    node_mock._run = mock._run
    node_mock._run.side_effect = exception("ARGH!")

    temp.mkdir()
    with pytest.raises(NodeError):
        node_mock.run(tmp_path)
    log_file = tmp_path / "xTMPx" / "pipe.errors"
    assert log_file.exists()

    error_text = log_file.read_text()
    assert "Errors =" in error_text
    assert paleomix.__version__ in error_text

    assert mock.mock_calls == [
        call._create_temp_dir(tmp_path),
        call._run(str(temp)),
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
    Node(**kwargs)._setup(None)


_INPUT_FILES_MISSING = (
    {"executables": ("ls", "shxxxx")},
    {"input_files": _OUT_FILES},
    {"auxiliary_files": _OUT_FILES},
)


@pytest.mark.parametrize("kwargs", _INPUT_FILES_MISSING)
def test__setup__input_files_missing(kwargs):
    with pytest.raises(NodeError):
        Node(**kwargs)._setup(None)


def test__teardown__output_files():
    Node(input_files=_EMPTY_FILE, output_files=_IN_FILES)._teardown(None)


def test__teardown__output_files_missing():
    node = Node(input_files=_EMPTY_FILE, output_files=_OUT_FILES)
    with pytest.raises(NodeError):
        node._teardown(None)


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


def test_node_collect_files__empty_folder(tmp_path: Path):
    assert list(Node._collect_files(tmp_path)) == []


def test_node_collect_files__root_files(tmp_path: Path):
    (tmp_path / "foo.txt").touch()
    (tmp_path / "bar.txt").touch()

    assert sorted(Node._collect_files(tmp_path)) == [
        "bar.txt",
        "foo.txt",
    ]


def test_node_collect_files__files_and_folders(tmp_path: Path):
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
    mock = _build_cmd_mock()

    node_mock = CommandNode(mock)
    node_mock._create_temp_dir = mock._test_node_._create_temp_dir
    node_mock._create_temp_dir.return_value = _DUMMY_TEMP
    node_mock._setup = mock._test_node_._setup
    node_mock._teardown = mock._test_node_._teardown
    node_mock._remove_temp_dir = mock._test_node_._remove_temp_dir

    node_mock.run(_DUMMY_TEMP_ROOT)

    assert mock.mock_calls == [
        call._test_node_._create_temp_dir(_DUMMY_TEMP_ROOT),
        call._test_node_._setup(_DUMMY_TEMP),
        call.run(_DUMMY_TEMP),
        call.join(),
        call._test_node_._teardown(_DUMMY_TEMP),
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


@pytest.mark.parametrize("kwargs", _SETUP_FILES_EXIST)
def test_commandnode_setup__files_exist(kwargs):
    cmd_mock = _build_cmd_mock(**kwargs)
    node = CommandNode(cmd_mock)
    node._setup(None)


_SETUP_FILES_MISSING = (
    {"executables": ("ls", "shxxxxxxxxxxx")},
    {"input_files": _OUT_FILES},
    {"auxiliary_files": _OUT_FILES},
)


@pytest.mark.parametrize("kwargs", _SETUP_FILES_MISSING)
def test_commandnode_setup__files_missing(kwargs):
    cmd_mock = _build_cmd_mock(**kwargs)
    node = CommandNode(cmd_mock)
    with pytest.raises(NodeError):
        node._setup(None)


###############################################################################
###############################################################################
# CommandNode: _run()


def test_commandnode_run__call_order():
    cmd_mock = _build_cmd_mock()
    node = CommandNode(cmd_mock)
    node._run("xTMPx")

    assert cmd_mock.mock_calls == [call.run("xTMPx"), call.join()]


def test_commandnode_run__exception_on_error():
    cmd_mock = _build_cmd_mock(return_codes=(1,))
    node = CommandNode(cmd_mock)
    with pytest.raises(CmdNodeError):
        node._run("xTMPx")

    assert cmd_mock.mock_calls == [call.run("xTMPx"), call.join()]


###############################################################################
###############################################################################
# CommandNode: _teardown


def _setup_temp_folders(tmp_path: Path):
    destination = tmp_path / "dst"
    tmp_path = tmp_path / "tmp"
    tmp_path.mkdir(parents=True, exist_ok=True)
    destination.mkdir(parents=True, exist_ok=True)
    return destination, tmp_path


# Commit is called on the command obj
def test_commandnode_teardown__commit(tmp_path: Path):
    cmd_mock = _build_cmd_mock()
    node = CommandNode(cmd_mock)
    node._teardown(tmp_path)
    assert cmd_mock.mock_calls == [call.commit(tmp_path)]


# Files exist in temp folder, and in destination after commit
def test_commandnode_teardown(tmp_path: Path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    cmd = AtomicCmd(
        ("echo", "-n", "1 2 3"),
        extra_files=[InputFile(_EMPTY_FILE)],
        stdout=str(destination / "foo.txt"),
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    assert (tmp_path / "foo.txt").exists()
    assert not (destination / "foo.txt").exists()
    node._teardown(tmp_path)
    assert not (tmp_path / "foo.txt").exists()
    assert (destination / "foo.txt").exists()


# Not all required files have been generated (atomic)
def test_commandnode_teardown__missing_files_in_temp(tmp_path: Path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    cmd = AtomicCmd(
        ("echo", "-n", "1 2 3"),
        extra_files=[
            InputFile(_EMPTY_FILE),
            OutputFile(str(destination / "bar.txt")),
        ],
        stdout=str(destination / "foo.txt"),
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    temp_files_before = set(tmp_path.iterdir())
    dest_files_before = set(destination.iterdir())

    with pytest.raises(CmdNodeError):
        node._teardown(tmp_path)
    assert temp_files_before == set(tmp_path.iterdir())
    assert dest_files_before == set(destination.iterdir())


# Not all specified TEMP_ files exist at _teardown (allowed)
def test_commandnode_teardown__missing_optional_files(tmp_path: Path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    cmd = AtomicCmd(
        ("echo", "-n", "1 2 3"),
        extra_files=[
            InputFile(_EMPTY_FILE),
            TempOutputFile("bar.txt"),
        ],
        stdout=str(destination / "foo.txt"),
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    node._teardown(tmp_path)
    assert os.listdir(str(tmp_path)) == []
    assert os.listdir(str(destination)) == ["foo.txt"]


# Not all required files were in place after commit
def test_commandnode_teardown__missing_files_in_dest(tmp_path: Path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    class _CmdMock(AtomicCmd):
        def commit(self, temp):
            AtomicCmd.commit(self, temp)
            (destination / "foo.txt").unlink()

    cmd = _CmdMock(
        (
            "touch",
            OutputFile(str(destination / "foo.txt")),
            OutputFile(str(destination / "bar.txt")),
        ),
        extra_files=[InputFile(_EMPTY_FILE)],
    )
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    node = CommandNode(cmd)
    with pytest.raises(NodeError):
        node._teardown(tmp_path)


def test_commandnode_teardown__extra_files_in_temp(tmp_path: Path):
    destination, tmp_path = _setup_temp_folders(tmp_path)

    unexpected_file = tmp_path / "unexpected_file.txt"
    unexpected_file.write_text("1 2 3")

    cmd = AtomicCmd(
        ("echo", "1 2 3"),
        extra_files=[InputFile(_EMPTY_FILE)],
        stdout=str(destination / "foo.txt"),
    )

    node = CommandNode(cmd)
    node._run(tmp_path)
    node._teardown(tmp_path)

    assert list(tmp_path.iterdir()) == [unexpected_file]
    assert list(destination.iterdir()) == [destination / "foo.txt"]
