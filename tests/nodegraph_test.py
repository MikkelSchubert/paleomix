# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

from paleomix.node import Node as Task
from paleomix.nodegraph import FileStatusCache, NodeGraph, StatusEnum


class MockCache(FileStatusCache):
    def __init__(
        self,
        abspaths: dict[str, str] | None = None,
        mtimes: dict[str, int] | None = None,
    ) -> None:
        self._abspaths = {} if abspaths is None else dict(abspaths)
        self._mtimes = {} if mtimes is None else dict(mtimes)

    def abspath(self, fpath: str) -> str:
        return self._abspaths.get(fpath, fpath)

    def _mtime_ns(self, fpath: str) -> int | None:
        return self._mtimes.get(fpath)

    def update_mtime(self, fpath: str, mtime: int) -> None:
        self._mtimes[fpath] = mtime


def single_task_status(task: Task) -> StatusEnum:
    return NodeGraph([task]).get_node_state(task)


def state_counts(**kwargs: int) -> dict[StatusEnum, int]:
    states = dict.fromkeys(StatusEnum, 0)
    for key, value in kwargs.items():
        states[StatusEnum(key)] += value

    return states


###############################################################################
###############################################################################


def test_empty_nodegraph() -> None:
    graph = NodeGraph(tasks=())
    assert graph.tasks == frozenset()
    assert graph.requirements == ()
    assert graph.get_state_counts() == dict.fromkeys(StatusEnum, 0)
    assert graph.get_and_reset_intermediate_files() == []
    assert graph.check_file_dependencies(MockCache())


def test_minimal_nodegraph_that_is_done() -> None:
    fscache = MockCache(mtimes={"input": 0, "output": 1})
    task = Task(input_files=["input"], output_files=["output"])
    graph = NodeGraph(tasks=[task], fscache=fscache)

    assert graph.tasks == frozenset([task])
    assert graph.check_file_dependencies(fscache)

    states = dict.fromkeys(StatusEnum, 0)
    states[StatusEnum.DONE] = 1
    assert graph.get_state_counts() == states


def test_minimal_nodegraph_that_is_ready() -> None:
    fscache = MockCache(mtimes={"input": 0})
    task = Task(input_files=["input"], output_files=["output"])
    graph = NodeGraph(tasks=[task], fscache=fscache)

    assert graph.tasks == frozenset([task])
    assert graph.get_state_counts() == state_counts(runnable=1)
    assert graph.check_file_dependencies(fscache)
    assert graph.get_state_counts() == state_counts(runnable=1)


def test_minimal_nodegraph_that_depends_on_missing_input_file() -> None:
    fscache = MockCache()
    task = Task(input_files=["input"], output_files=["output"])
    graph = NodeGraph(tasks=[task], fscache=fscache)

    assert graph.get_node_state(task) == StatusEnum.RUNNABLE
    assert graph.get_state_counts() == state_counts(runnable=1)
    assert not graph.check_file_dependencies(fscache)
    assert graph.get_node_state(task) == StatusEnum.ERROR
    assert graph.get_state_counts() == state_counts(error=1)

    states = dict.fromkeys(StatusEnum, 0)
    states[StatusEnum.ERROR] = 1
    assert graph.get_state_counts() == states


def test_minimal_nodegraph_that_depends_on_missing_auxiliary_file() -> None:
    fscache = MockCache(mtimes={"input": 0})
    task = Task(input_files=["input"], auxiliary_files=["aux"], output_files=["output"])
    graph = NodeGraph(tasks=[task], fscache=fscache)

    assert graph.get_node_state(task) == StatusEnum.RUNNABLE
    assert graph.get_state_counts() == state_counts(runnable=1)
    assert not graph.check_file_dependencies(fscache)
    assert graph.get_node_state(task) == StatusEnum.ERROR
    assert graph.get_state_counts() == state_counts(error=1)


def test_check_empty_version_requirements_is_ok() -> None:
    assert NodeGraph.check_version_requirements(())
    assert NodeGraph.check_version_requirements((), force=True)


# TODO: Detect tasks that only generate intermediate files, none of which are used
#       by a downstream task, as these will not be run (unless strategy = required)
