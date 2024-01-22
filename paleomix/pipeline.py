#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import logging
import math
import multiprocessing
import os
import signal
import sys
import time
from shlex import quote
from typing import TYPE_CHECKING, Any, Iterable, NoReturn

import paleomix.common.logging
import paleomix.core.reports
from paleomix.common.fileutils import try_remove
from paleomix.common.procs import terminate_all_processes
from paleomix.common.text import format_timespan, padded_table
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.core.workers import EVT_CAPACITY, EVT_SHUTDOWN, EVT_TASK_DONE, Manager
from paleomix.node import Node, NodeError, NodeMissingFilesError
from paleomix.nodegraph import (
    CleanupStrategy,
    FileStatusCache,
    NodeGraph,
    NodeGraphError,
    StatusEnum,
)

if TYPE_CHECKING:
    from configargparse import ArgumentParser

    from paleomix.common.argparse import ArgumentGroup


class Pypeline:
    def __init__(
        self,
        nodes: Iterable[Node],
        temp_root: str = "/tmp",
        max_threads: int = 1,
        intermediate_files: CleanupStrategy = CleanupStrategy.DELETE,
        required_files: Iterable[str] = (),
    ) -> None:
        self._nodes = safe_coerce_to_tuple(nodes)
        for node in self._nodes:
            if not isinstance(node, Node):
                raise TypeError("Node object expected, received %s" % repr(node))

        self._logger = logging.getLogger(__name__)
        # Set if a keyboard-interrupt (SIGINT) has been caught
        self._interrupted = 0
        self._threads = max(0, max_threads)
        self._temp_root = temp_root
        self._intermediate_files_strategy = intermediate_files
        self._required_files = frozenset(required_files)
        self._start_times: dict[Node, float] = {}
        self._progress_color: str | None = None

        self._event_handlers = {
            # The number of available threads has changed
            EVT_CAPACITY: self._event_capacity,
            # A task has finished running
            EVT_TASK_DONE: self._event_task_done,
            # Worker was shut down, possibly killing tasks in the process
            EVT_SHUTDOWN: self._event_shutdown,
        }

    def run(self, mode: str = "run") -> int:
        fscache = FileStatusCache()

        try:
            nodegraph = NodeGraph(
                tasks=self._nodes,
                fscache=fscache,
                intermediate_files=self._intermediate_files_strategy,
                required_files=self._required_files,
            )
        except NodeGraphError as error:
            self._logger.error(error)
            return 1

        if mode not in ("run", "dry_run"):
            return self._print_report(mode, nodegraph, fscache)
        elif not nodegraph.check_file_dependencies(fscache):
            return 1

        manager = Manager(
            threads=self._threads,
            temp_root=self._temp_root,
            requirements=nodegraph.requirements,
        )

        try:
            # Handle setup/teardown of command-line interface and termination of workers
            with manager:
                if not manager.start():
                    self._logger.error("Manager failed to start; terminating")
                    return 1

                if mode == "dry_run":
                    # Wait for any remote workers discovered during startup
                    if not manager.wait_for_workers():
                        self._logger.error("Workers failed to start; terminating")
                        return 1

                    self._summarize_pipeline(nodegraph, dry_run=mode == "dry_run")
                    return 0

                self._logger.info("Running pipeline (press 'h' for help):")

                # Install signal handlers to allow graceful termination
                signal.signal(signal.SIGINT, self._sigint_handler)
                signal.signal(signal.SIGHUP, self._sigterm_handler)
                signal.signal(signal.SIGTERM, self._sigterm_handler)

                return self._run(nodegraph, manager)
        finally:
            terminate_all_processes()

            for filename in paleomix.common.logging.get_logfiles():
                self._logger.info("Log-file written to %r", filename)

    def _run(self, nodegraph: NodeGraph, manager: Manager) -> int:
        # Set of remaining nodes to be run
        tasks: dict[Node, _TaskInfo] = {}
        for task in nodegraph.tasks:
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

            self._clean_intermediate_files(nodegraph)

        self._clean_intermediate_files(nodegraph)

        self._logger.info("Shutting down workers")
        manager.shutdown()

        self._summarize_pipeline(nodegraph, verbose=any_errors)

        return 1 if any_errors else 0

    def _event_capacity(
        self,
        nodegraph: NodeGraph,
        manager: Manager,
        tasks: dict[Node, _TaskInfo],
        worker: str,
        **event: Any,
    ) -> bool:
        if not self._interrupted:
            idle_threads = event["threads"]
            for task, task_info in sorted(tasks.items(), key=lambda it: it[0].id):
                if worker in task_info.blacklisted_from:
                    continue

                if nodegraph.get_node_state(task) == nodegraph.RUNNABLE:
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
        tasks: dict[Node, _TaskInfo],
        worker: str,
        **event: Any,
    ) -> bool:
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
                self._set_node_state(nodegraph, task, nodegraph.RUNNABLE)
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
        tasks: dict[Node, _TaskInfo],
        worker: str,
        worker_name: str,
        **event: Any,
    ):
        self._logger.error("PALEOMIX worker %s terminated", worker_name)

        any_errors = False
        workers = manager.workers.keys()
        for task, task_info in tuple(tasks.items()):
            if task_info.running_on == worker:
                self._logger.warning("Re-trying %s", task)
                self._set_node_state(nodegraph, task, nodegraph.RUNNABLE)
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

    def _prune_tasks(self, nodegraph: NodeGraph, tasks: dict[Node, _TaskInfo]):
        # The completion or failure of a task may result in the failure/completion of
        # any number of other tasks, the latter when tasks depend on validation steps
        for task in tuple(tasks):
            if nodegraph.get_node_state(task) in (nodegraph.DONE, nodegraph.ERROR):
                tasks.pop(task)

    def _clean_intermediate_files(self, nodegraph: NodeGraph) -> None:
        # FIXME: Should not block main thread
        for filepath in nodegraph.get_and_reset_intermediate_files():
            self._logger.debug("Removing no longer needed file %s", quote(filepath))
            if not try_remove(filepath):
                self._logger.error("Failed to remove temp file %s", quote(filepath))

    def _handle_task_error(
        self,
        nodegraph: NodeGraph,
        task: Node,
        error: Any,
        backtrace: list[str] | None,
        **kwargs: Any,
    ):
        self._set_node_state(nodegraph, task, nodegraph.ERROR)

        if not isinstance(error, NodeError):
            error = f"Unhandled exception while running {task}:"

        message = str(error).split("\n")

        if backtrace:
            message.extend("".join(backtrace).rstrip().split("\n"))

        if isinstance(error, NodeError) and error.path:
            message.append("For more information about this error, see")
            message.append("  " + quote(os.path.join(error.path, "pipe.errors")))

        self._logger.error("\n".join(message))

    def _print_report(
        self,
        mode: str,
        graph: NodeGraph,
        fscache: FileStatusCache,
    ) -> int:
        try:
            if mode == "input_files":
                self._logger.info("Collecting and printing input files ..")
                return paleomix.core.reports.input_files(graph, fscache)
            elif mode == "output_files":
                self._logger.info("Collecting and printing output files ..")
                return paleomix.core.reports.output_files(graph, fscache)
            elif mode == "executables":
                self._logger.info("Collecting and printing required executables ..")
                return paleomix.core.reports.required_executables(graph)
            elif mode == "pipeline_tasks":
                self._logger.info("Printing pipeline tasks ..")
                return paleomix.core.reports.pipeline_tasks(graph)
            else:
                raise ValueError(f"Unknown pipeline mode {mode!r}")
        except BrokenPipeError:
            return 0
        except NodeGraphError as error:
            self._logger.error(error)
            return 1

    def _sigint_handler(self, signum: int, frame: object) -> None:
        """Signal handler; see signal.signal."""
        now = time.time()
        if not self._interrupted:
            self._interrupted = now
            self._logger.warning(
                "Keyboard interrupt detected, waiting for running tasks to complete. "
                "Press CTRL-C again within the next 5 seconds to terminate immediately."
            )
        elif now - self._interrupted <= 5.0:
            self._sigterm_handler(signum, frame)
        else:
            self._interrupted = now
            self._logger.warning(
                "Pipeline is waiting for running tasks to terminate. Press CTRL-C "
                "again within the next 5 seconds to terminate immediately."
            )

    def _sigterm_handler(self, signum: int, frame: object) -> NoReturn:
        self._logger.warning("Terminating due to signal %i", signum)

        signal.signal(signal.SIGINT, signal.SIG_DFL)
        signal.signal(signal.SIGHUP, signal.SIG_DFL)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)

        sys.exit(-signum)

    def _summarize_pipeline(
        self,
        nodegraph: NodeGraph,
        *,
        dry_run: bool = False,
        verbose: bool = True,
    ) -> None:
        states = nodegraph.get_state_counts()

        if verbose:
            rows = [
                ("Number of tasks:", sum(states.values())),
                ("Number of done tasks:", states[nodegraph.DONE]),
                ("Number of runnable tasks:", states[nodegraph.RUNNABLE]),
                ("Number of queued tasks:", states[nodegraph.QUEUED]),
                ("Number of failed tasks:", states[nodegraph.ERROR]),
            ]

            for message in padded_table(rows):
                self._logger.info(message)

        if states[nodegraph.ERROR]:
            self._logger.warning("Errors were detected in pipeline")
        elif self._interrupted:
            self._logger.info("Pipeline interrupted by user")
        elif dry_run:
            self._logger.info("Dry run completed successfully")
        else:
            self._logger.info("Pipeline completed successfully")

    def _set_node_state(self, graph: NodeGraph, node: Node, state: StatusEnum) -> None:
        for node, old_state, new_state in graph.set_node_status(node, state):
            if new_state in (NodeGraph.RUNNING, NodeGraph.DONE):
                runtime = ""
                if new_state == NodeGraph.RUNNING:
                    self._start_times[node] = time.time()
                    event = "Started"
                elif old_state == NodeGraph.RUNNING:
                    event = "Finished"
                    end_time = time.time()
                    start_time = self._start_times.pop(node)
                    runtime = f" in {format_timespan(end_time - start_time)}"
                else:
                    event = "Already finished"

                status = _Progress(graph.get_state_counts(), self._progress_color)
                extra = {"status": status}

                self._logger.info("%s %s%s", event, node, runtime, extra=extra)
            elif new_state == graph.ERROR:
                self._progress_color = "red"


def add_argument_groups(parser: ArgumentParser) -> None:
    add_scheduling_argument_group(parser)
    add_io_argument_group(parser)


def add_scheduling_argument_group(parser: ArgumentParser) -> ArgumentGroup:
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

    group = parser.add_argument_group("Intermediate files")
    group.add_argument(
        "--intermediate-files",
        type=CleanupStrategy,
        default=CleanupStrategy.DELETE,
        choices=tuple(CleanupStrategy),
        help="Strategy for handling intermediate files: 'delete' removes intermediate "
        "as soon as they are no longer needed; 'keep' keeps any existing intermediate "
        "files, but does not regenerate missing intermediate files; and 'require' will "
        "re-run tasks if intermediate files are missing",
    )
    group.add_argument(
        "--require-files",
        metavar="glob",
        default=[],
        action="append",
        help="Intermediate files matching the glob will never be deleted and will be "
        "re-generated if missing, no matter what --intermediate-files strategy is used",
    )

    return group


def add_io_argument_group(parser: ArgumentParser) -> ArgumentGroup:
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
    def __init__(self, task: Node) -> None:
        self.task = task
        self.running_on: str | None = None
        self.blacklisted_from: dict[str, object] = {}
        self.start_time = 0.0


class _Progress(paleomix.common.logging.Status):
    def __init__(self, state_counts: dict[StatusEnum, int], color: str | None) -> None:
        super().__init__(color)
        self._state_counts = state_counts

    def __str__(self) -> str:
        total = sum(self._state_counts.values())
        nth = self._state_counts[NodeGraph.DONE] + self._state_counts[NodeGraph.ERROR]

        if total > 200:
            value = f"{math.floor((1000 * nth) / total) / 10: >5.1f}%"
        elif total > 0:
            value = f"{(100 * nth) // total: >3g}%"
        else:
            value = "N/A"

        return value
