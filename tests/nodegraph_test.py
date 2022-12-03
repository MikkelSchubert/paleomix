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
from typing import Dict, Optional

from paleomix.node import Node as Task
from paleomix.nodegraph import FileStatusCache, NodeGraph, StatusEnum


class MockCache(FileStatusCache):
    def __init__(
        self,
        abspaths: Optional[Dict[str, str]] = None,
        mtimes: Optional[Dict[str, int]] = None,
    ):
        self._abspaths = {} if abspaths is None else dict(abspaths)
        self._mtimes = {} if mtimes is None else dict(mtimes)

    def abspath(self, fpath: str) -> str:
        return self._abspaths.get(fpath, fpath)

    def _mtime_ns(self, fpath: str) -> Optional[int]:
        return self._mtimes.get(fpath)

    def update_mtime(self, fpath: str, mtime: int) -> None:
        self._mtimes[fpath] = mtime


def single_task_status(task: Task) -> StatusEnum:
    return NodeGraph([task]).get_node_state(task)


def state_counts(**kwargs: int) -> Dict[StatusEnum, int]:
    states = dict.fromkeys(StatusEnum, 0)
    for key, value in kwargs.items():
        states[StatusEnum(key)] += value

    return states


###############################################################################
###############################################################################


def test_empty_nodegraph():
    graph = NodeGraph(tasks=())
    assert graph.tasks == frozenset()
    assert graph.requirements == ()
    assert graph.get_state_counts() == dict.fromkeys(StatusEnum, 0)
    assert graph.get_and_reset_intermediate_files() == []
    assert graph.check_file_dependencies(MockCache())


def test_minimal_nodegraph_that_is_done():
    fscache = MockCache(mtimes={"input": 0, "output": 1})
    task = Task(input_files=["input"], output_files=["output"])
    graph = NodeGraph(tasks=[task], fscache=fscache)

    assert graph.tasks == frozenset([task])
    assert graph.check_file_dependencies(fscache)

    states = dict.fromkeys(StatusEnum, 0)
    states[StatusEnum.DONE] = 1
    assert graph.get_state_counts() == states


def test_minimal_nodegraph_that_is_ready():
    fscache = MockCache(mtimes={"input": 0})
    task = Task(input_files=["input"], output_files=["output"])
    graph = NodeGraph(tasks=[task], fscache=fscache)

    assert graph.tasks == frozenset([task])
    assert graph.get_state_counts() == state_counts(runnable=1)
    assert graph.check_file_dependencies(fscache)
    assert graph.get_state_counts() == state_counts(runnable=1)


def test_minimal_nodegraph_that_depends_on_missing_input_file():
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


def test_minimal_nodegraph_that_depends_on_missing_auxiliary_file():
    fscache = MockCache(mtimes={"input": 0})
    task = Task(input_files=["input"], auxiliary_files=["aux"], output_files=["output"])
    graph = NodeGraph(tasks=[task], fscache=fscache)

    assert graph.get_node_state(task) == StatusEnum.RUNNABLE
    assert graph.get_state_counts() == state_counts(runnable=1)
    assert not graph.check_file_dependencies(fscache)
    assert graph.get_node_state(task) == StatusEnum.ERROR
    assert graph.get_state_counts() == state_counts(error=1)


def test_check_empty_version_requirements_is_ok():
    assert NodeGraph.check_version_requirements(())
    assert NodeGraph.check_version_requirements((), force=True)


# TODO: Detect tasks that only generate intermediate files, none of which are used
#       by a downstream task, as these will not be run (unless strategy = required)
