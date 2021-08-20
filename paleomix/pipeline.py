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
import argparse
import logging
import math
import multiprocessing
import os
import signal
import sys
import time
from shlex import quote
from typing import IO, Any, Dict, Iterable, List, Optional

import paleomix.common.logging
import paleomix.core.reports
from paleomix.common.text import format_timespan, padded_table
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.core.workers import EVT_CAPACITY, EVT_SHUTDOWN, EVT_TASK_DONE, Manager
from paleomix.node import Node, NodeError, NodeMissingFilesError
from paleomix.nodegraph import NodeGraph, NodeGraphError


class Pypeline:
    def __init__(
        self,
        nodes: Iterable[Node],
        temp_root: str = "/tmp",
        max_threads: int = 1,
        implicit_dependencies: bool = False,
    ):
        self._nodes = safe_coerce_to_tuple(nodes)
        for node in self._nodes:
            if not isinstance(node, Node):
                raise TypeError("Node object expected, recieved %s" % repr(node))

        self._logger = logging.getLogger(__name__)
        # Set if a keyboard-interrupt (SIGINT) has been caught
        self._interrupted = False
        self._threads = max(0, max_threads)
        self._temp_root = temp_root
        self._implicit_dependencies = implicit_dependencies
        self._start_times: Dict[Node, float] = {}
        self._progress_color: Optional[str] = None

        self._event_handlers = {
            # The number of available threads has changed
            EVT_CAPACITY: self._event_capacity,
            # A task has finished running
            EVT_TASK_DONE: self._event_task_done,
            # Worker was shut down, possibily killing tasks in the process
            EVT_SHUTDOWN: self._event_shutdown,
        }

    def run(self, mode: str = "run") -> int:
        if mode == "dry_run":
            self._logger.info("Dry running pipeline ..")
        elif mode == "run":
            self._logger.info("Starting pipeline ..")
        else:
            return self._print_report(mode)

        try:
            nodegraph = NodeGraph(
                nodes=self._nodes,
                implicit_dependencies=self._implicit_dependencies,
            )
        except NodeGraphError as error:
            self._logger.error(error)
            return 1

        manager = Manager(
            threads=self._threads,
            temp_root=self._temp_root,
            requirements=nodegraph.requirements,
        )

        try:
            # Handle setup/teardown of commandline interface and termination of workers
            with manager:
                if not manager.start():
                    self._logger.error("Manager failed to start; terminating")
                    return 1

                if mode == "dry_run":
                    # Wait for any remote workers discovered during startup
                    if not manager.wait_for_workers():
                        self._logger.error("Workers failed to start; terminating")
                        return 1

                    self._summarize_pipeline(nodegraph)
                    self._logger.info("Dry run done")
                    return 0

                # Install signal handlers to allow graceful termination
                signal.signal(signal.SIGINT, self._sigint_handler)
                signal.signal(signal.SIGHUP, self._sigterm_handler)
                signal.signal(signal.SIGTERM, self._sigterm_handler)

                return self._run(nodegraph, manager)
        finally:
            for filename in paleomix.common.logging.get_logfiles():
                self._logger.info("Log-file written to %r", filename)

    def _run(self, nodegraph: NodeGraph, manager: Manager) -> int:
        # Set of remaining nodes to be run
        tasks: Dict[Node, _TaskInfo] = {}
        for task in nodegraph.iterflat():
            state = nodegraph.get_node_state(task)
            if state not in (nodegraph.DONE, nodegraph.ERROR):
                tasks[task] = _TaskInfo(task)

        any_errors = False
        # Keep looping as long as there are tasks left to run or tasks running
        while (tasks and not self._interrupted) or any(manager.tasks):
            for event in manager.poll():
                handler = self._event_handlers.get(event["event"])
                if handler is None:
                    self._logger.error("Unknown event in pipeline: %r", event)
                elif not handler(nodegraph, manager, tasks, **event):
                    any_errors = True

        self._logger.info("Shutting down workers")
        manager.shutdown()

        self._summarize_pipeline(nodegraph, verbose=any_errors)

        return 1 if any_errors else 0

    def _event_capacity(
        self,
        nodegraph: NodeGraph,
        manager: Manager,
        tasks: Dict[Node, "_TaskInfo"],
        worker: str,
        **event: Any
    ):
        if not self._interrupted:
            idle_threads = event["threads"]
            for task, task_info in sorted(tasks.items(), key=lambda it: it[0].id):
                if worker in task_info.blacklisted_from:
                    continue

                if nodegraph.get_node_state(task) == nodegraph.RUNABLE:
                    if idle_threads >= task.threads or event["overcommit"]:
                        if not manager.start_task(worker, task):
                            # Error in worker; this will probably be picked up next loop
                            return False

                        task_info.running_on = worker
                        self._set_node_state(nodegraph, task, nodegraph.RUNNING)

                        # Overcommiting allowed only if a worker is idle
                        event["overcommit"] = False

                        idle_threads -= task.threads
                        if idle_threads <= 0:
                            break

        return True

    def _event_task_done(
        self,
        nodegraph: NodeGraph,
        manager: Manager,
        tasks: Dict[Node, "_TaskInfo"],
        worker: str,
        **event: Any
    ):
        any_errors = False
        task = event["task"]
        task_info = tasks.pop(task)
        if event["error"] is None:
            self._set_node_state(nodegraph, task, nodegraph.DONE)
        elif isinstance(event["error"], NodeMissingFilesError):
            # Node was unexpectedly missing input files; this is either a programming
            # error, a user deleting stuff, or NFS (caches) not having been updated, so
            # we try to run it on another worker, if any are available.
            task_info.blacklisted_from[worker] = event
            if manager.workers.keys() - task_info.blacklisted_from:
                self._logger.warning("Re-trying %s", task)
                self._set_node_state(nodegraph, task, nodegraph.RUNABLE)
                tasks[task] = task_info
            else:  # No more nodes left to try so we'll just have to error out
                self._handle_task_error(nodegraph, **event)
                any_errors = True
        else:  # Permanent failure
            self._handle_task_error(nodegraph, **event)
            any_errors = True

        self._prune_tasks(nodegraph, tasks)

        return not any_errors

    def _event_shutdown(
        self,
        nodegraph: NodeGraph,
        manager: Manager,
        tasks: Dict[Node, "_TaskInfo"],
        worker: str,
        worker_name: str,
        **event: Any
    ):
        self._logger.error("PALEOMIX worker %s terminated", worker_name)

        any_errors = False
        workers = manager.workers.keys()
        for task, task_info in tuple(tasks.items()):
            if task_info.running_on == worker:
                self._logger.warning("Re-trying %s", task)
                self._set_node_state(nodegraph, task, nodegraph.RUNABLE)
                task_info.running_on = None

            # Check nodes that can no longer be completed
            if not (workers - task_info.blacklisted_from):
                # Pick arbitrary error message
                for event in task_info.blacklisted_from.values():
                    self._handle_task_error(nodegraph, **event)
                    any_errors = True
                    tasks.pop(task)

        if any_errors:
            self._prune_tasks(nodegraph, tasks)

        return not any_errors

    def _prune_tasks(self, nodegraph: NodeGraph, tasks: Dict[Node, "_TaskInfo"]):
        # The completion or failure of a task may result in the failure/completion of
        # any number of other tasks, the latter when tasks depend on validation steps
        for task in tuple(tasks):
            if nodegraph.get_node_state(task) in (nodegraph.DONE, nodegraph.ERROR):
                tasks.pop(task)

    def _handle_task_error(
        self,
        nodegraph: NodeGraph,
        task: Node,
        error: Any,
        backtrace: Optional[List[str]],
        **kwargs: Any
    ):
        self._set_node_state(nodegraph, task, nodegraph.ERROR)

        if not isinstance(error, NodeError):
            error = "Unhandled exception while running {}:".format(task)

        message = str(error).split("\n")

        if backtrace:
            message.extend("".join(backtrace).rstrip().split("\n"))

        if isinstance(error, NodeError) and error.path:
            message.append("For more information about this error, see")
            message.append("  " + quote(os.path.join(error.path, "pipe.errors")))

        self._logger.error("\n".join(message))

    def _print_report(self, mode: str, file: IO[str] = sys.stdout) -> int:
        try:
            if mode == "input_files":
                self._logger.info("Collecting and printing input files ..")
                return paleomix.core.reports.input_files(self._nodes, file)
            elif mode == "output_files":
                self._logger.info("Collecting and printing output files ..")
                return paleomix.core.reports.output_files(self._nodes, file)
            elif mode == "executables":
                self._logger.info("Collecting and printing required executables ..")
                return paleomix.core.reports.required_executables(self._nodes, file)
            elif mode == "pipeline_tasks":
                self._logger.info("Printing pipeline tasks ..")
                return paleomix.core.reports.pipeline_tasks(self._nodes, file)
            else:
                raise ValueError("Unknown pipeline mode {!r}".format(mode))
        except BrokenPipeError:
            return 0
        except NodeGraphError as error:
            self._logger.error(error)
            return 1

    def _sigint_handler(self, signum: int, frame: Any):
        """Signal handler; see signal.signal."""
        if not self._interrupted:
            self._interrupted = True
            self._logger.warning(
                "Keyboard interrupt detected, waiting for running tasks to complete. "
                "Press CTRL-C again to force termination."
            )
        else:
            self._sigterm_handler(signum, frame)

    def _sigterm_handler(self, signum: int, frame: Any):
        self._logger.warning("Terminating due to signal %i", signum)
        sys.exit(-signum)

    def _summarize_pipeline(self, nodegraph: NodeGraph, verbose: bool = True):
        states = nodegraph.get_state_counts()

        if verbose:
            rows = [
                ("Number of tasks:", sum(states)),
                ("Number of done tasks:", states[nodegraph.DONE]),
                ("Number of runable tasks:", states[nodegraph.RUNABLE]),
                ("Number of queued tasks:", states[nodegraph.QUEUED]),
                ("Number of outdated tasks:", states[nodegraph.OUTDATED]),
                ("Number of failed tasks:", states[nodegraph.ERROR]),
            ]

            for message in padded_table(rows):
                self._logger.info(message)

        if states[nodegraph.ERROR]:
            self._logger.warning("Errors were detected while running pipeline")
        else:
            self._logger.info("Pipeline completed successfully")

    def _set_node_state(self, graph: NodeGraph, node: Node, state: int) -> None:
        for node, old_state, new_state in graph.set_node_state(node, state):
            if new_state in (NodeGraph.RUNNING, NodeGraph.DONE):
                runtime = ""
                if new_state == NodeGraph.RUNNING:
                    self._start_times[node] = time.time()
                    event = "Started"
                elif old_state == NodeGraph.RUNNING:
                    event = "Finished"
                    end_time = time.time()
                    start_time = self._start_times.pop(node)
                    runtime = " in {}".format(format_timespan(end_time - start_time))
                else:
                    event = "Already finished"

                status = _Progress(graph.get_state_counts(), self._progress_color)
                extra = {"status": status}

                self._logger.info("%s %s%s", event, node, runtime, extra=extra)
            elif new_state == graph.ERROR:
                self._progress_color = "red"


def add_argument_groups(parser: argparse.ArgumentParser) -> None:
    add_scheduling_argument_group(parser)
    add_io_argument_group(parser)


def add_scheduling_argument_group(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:
    parser.set_defaults(pipeline_mode="run")

    group = parser.add_argument_group("Pipeline Scheduling")
    group.add_argument(
        "--dry-run",
        const="dry_run",
        action="store_const",
        dest="pipeline_mode",
        help="Build pipeline and check prerequisites, but do not execute any tasks",
    )
    group.add_argument(
        "--max-threads",
        type=int,
        default=max(2, multiprocessing.cpu_count()),
        help="Max number of threads to use in total",
    )

    return group


def add_io_argument_group(parser: argparse.ArgumentParser) -> argparse._ArgumentGroup:
    parser.set_defaults(pipeline_mode="run")

    group = parser.add_argument_group("Pipeline Input/Output")
    group.add_argument(
        "--list-input-files",
        const="input_files",
        action="store_const",
        dest="pipeline_mode",
        help="List all external input files used by pipeline for the makefile(s)",
    )
    group.add_argument(
        "--list-output-files",
        const="output_files",
        action="store_const",
        dest="pipeline_mode",
        help="List all output files generated by pipeline for the makefile(s)",
    )
    group.add_argument(
        "--list-executables",
        const="executables",
        action="store_const",
        dest="pipeline_mode",
        help="List executables used by the pipeline and their version requirements.",
    )
    group.add_argument(
        "--list-pipeline-tasks",
        const="pipeline_tasks",
        action="store_const",
        dest="pipeline_mode",
        help="Print tree of tasks in the current pipeline",
    )

    return group


class _TaskInfo:
    def __init__(self, task: Node):
        self.task = task
        self.running_on: Optional[str] = None
        self.blacklisted_from: Dict[str, Any] = {}
        self.start_time = 0.0


class _Progress(paleomix.common.logging.Status):
    def __init__(self, state_counts: List[int], color: Optional[str]):
        super().__init__(color)
        self._state_counts = state_counts

    def __str__(self):
        total = sum(self._state_counts)
        nth = self._state_counts[NodeGraph.DONE] + self._state_counts[NodeGraph.ERROR]

        if total > 200:
            value = "{: >5.1f}%".format(math.floor((1000 * nth) / total) / 10)
        elif total > 0:
            value = "{: >3g}%".format((100 * nth) // total)
        else:
            value = "N/A"

        return value
