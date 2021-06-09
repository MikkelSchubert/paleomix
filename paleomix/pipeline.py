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
import logging
import multiprocessing
import os
import signal
import sys

from shlex import quote

import paleomix.common.logging
import paleomix.core.reports

from paleomix.common.text import padded_table
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.core.workers import (
    EVT_CAPACITY,
    EVT_SHUTDOWN,
    EVT_TASK_DONE,
    Manager,
)
from paleomix.node import Node, NodeError, NodeMissingFilesError
from paleomix.nodegraph import NodeGraph, NodeGraphError


class Pypeline:
    def __init__(
        self,
        nodes,
        temp_root="/tmp",
        max_threads=1,
        implicit_dependencies=False,
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
        self._manager = None

        self._event_handlers = {
            # The number of available threads has changed
            EVT_CAPACITY: self._event_capacity,
            # A task has finished running
            EVT_TASK_DONE: self._event_task_done,
            # Worker was shut down, possibily killing tasks in the process
            EVT_SHUTDOWN: self._event_shutdown,
        }

    def run(self, mode="run"):
        if mode not in ("run", "dry_run", "input_files", "output_files", "executables"):
            raise ValueError("Unknown pipeline mode {!r}".format(mode))
        elif mode == "dry_run":
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

        self._manager = Manager(
            threads=self._threads,
            temp_root=self._temp_root,
            requirements=nodegraph.requirements,
        )

        sigint_handler = signal.getsignal(signal.SIGINT)

        try:
            # Handle setup/teardown of commandline interface and termination of workers
            with self._manager:
                if not self._manager.start():
                    self._logger.error("Manager failed to start; terminating")
                    return 1

                if mode == "dry_run":
                    # Wait for any remote workers discovered during startup
                    if not self._manager.wait_for_workers():
                        self._logger.error("Workers failed to start; terminating")
                        return 1

                    self._summarize_pipeline(nodegraph)
                    self._logger.info("Dry run done")
                    return 0

                # Install signal handler used to allow graceful termination
                signal.signal(signal.SIGINT, self._sigint_handler)

                return self._run(nodegraph)
        finally:
            signal.signal(signal.SIGINT, sigint_handler)
            for filename in paleomix.common.logging.get_logfiles():
                self._logger.info("Log-file written to %r", filename)

    def _run(self, nodegraph):
        # Set of remaining nodes to be run
        tasks = {}
        for task in nodegraph.iterflat():
            state = nodegraph.get_node_state(task)
            if state not in (nodegraph.DONE, nodegraph.ERROR):
                tasks[task] = {
                    "running_on": None,
                    "blacklisted_from": {},
                }

        any_errors = False
        # Keep looping as long as there are tasks left to run or tasks running
        while (tasks and not self._interrupted) or any(self._manager.tasks):
            for event in self._manager.poll():
                handler = self._event_handlers.get(event["event"])
                if handler is None:
                    self._logger.error("Unknown event in pipeline: %r", event)
                elif not handler(nodegraph, tasks, **event):
                    any_errors = True

        self._logger.info("Shutting down workers")
        self._manager.shutdown()

        self._summarize_pipeline(nodegraph, verbose=any_errors)

        return 1 if any_errors else 0

    def _event_capacity(self, nodegraph, tasks, worker, **event):
        if not self._interrupted:
            idle_threads = event["threads"]
            for task, task_info in sorted(tasks.items(), key=lambda it: it[0].id):
                if worker in task_info["blacklisted_from"]:
                    continue

                if nodegraph.get_node_state(task) == nodegraph.RUNABLE:
                    if idle_threads >= task.threads or event["overcommit"]:
                        if not self._manager.start_task(worker, task):
                            # Error in worker; this will probably be picked up next loop
                            return False

                        task_info["running_on"] = worker
                        nodegraph.set_node_state(task, nodegraph.RUNNING)

                        # Overcommiting allowed only if a worker is idle
                        event["overcommit"] = False

                        idle_threads -= task.threads
                        if idle_threads <= 0:
                            break

        return True

    def _event_task_done(self, nodegraph, tasks, worker, **event):
        any_errors = False
        task = event["task"]
        task_info = tasks.pop(task)
        if event["error"] is None:
            nodegraph.set_node_state(task, nodegraph.DONE)
        elif isinstance(event["error"], NodeMissingFilesError):
            # Node was unexpectedly missing input files; this is either a programming
            # error, a user deleting stuff, or NFS (caches) not having been updated, so
            # we try to run it on another worker, if any are available.
            task_info["blacklisted_from"][worker] = event
            if self._manager.workers - task_info["blacklisted_from"]:
                self._logger.warning("Re-trying %s", task)
                nodegraph.set_node_state(task, nodegraph.RUNABLE)
                tasks[task] = task_info
            else:  # No more nodes left to try so we'll just have to error out
                self._handle_task_error(nodegraph, **event)
                any_errors = True
        else:  # Permanent failure
            self._handle_task_error(nodegraph, **event)
            any_errors = True

        self._prune_tasks(nodegraph, tasks)

        return not any_errors

    def _event_shutdown(self, nodegraph, tasks, worker, worker_name, **event):
        self._logger.error("PALEOMIX worker %s terminated", worker_name)

        any_errors = False
        workers = self._manager.workers
        for task, task_info in tuple(tasks.items()):
            if task_info["running_on"] == worker:
                self._logger.warning("Re-trying %s", task)
                nodegraph.set_node_state(task, nodegraph.RUNABLE)
                task_info["running_on"] = None

            # Check nodes that can no longer be completed
            if not (workers - task_info["blacklisted_from"]):
                # Pick arbitrary error message
                for event in task_info["blacklisted_from"].values():
                    self._handle_task_error(nodegraph, **event)
                    any_errors = True
                    tasks.pop(task)

        if any_errors:
            self._prune_tasks(nodegraph, tasks)

        return not any_errors

    def _prune_tasks(self, nodegraph, tasks):
        # The completion or failure of a task may result in the failure/completion of
        # any number of other tasks, the latter when tasks depend on validation steps
        for task in tuple(tasks):
            if nodegraph.get_node_state(task) in (nodegraph.DONE, nodegraph.ERROR):
                tasks.pop(task)

    def _handle_task_error(self, nodegraph, task, error, backtrace, **kwargs):
        nodegraph.set_node_state(task, nodegraph.ERROR)

        if not isinstance(error, NodeError):
            error = "Unhandled exception while running {}:".format(task)

        message = str(error).split("\n")

        if backtrace:
            message.extend("".join(backtrace).rstrip().split("\n"))

        if isinstance(error, NodeError) and error.path:
            message.append("For more information about this error, see")
            message.append("  " + quote(os.path.join(error.path, "pipe.errors")))

        self._logger.error("\n".join(message))

    def _print_report(self, mode, file=sys.stdout):
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
            else:
                raise ValueError(mode)
        except BrokenPipeError:
            return 0
        except NodeGraphError as error:
            self._logger.error(error)
            return 1

    def _sigint_handler(self, signum, frame):
        """Signal handler; see signal.signal."""
        if not self._interrupted:
            self._interrupted = True
            self._logger.warning(
                "Keyboard interrupt detected, waiting for running tasks to complete. "
                "Press CTRL-C again to force termination."
            )
        else:
            self._logger.warning("Terminating pipeline!")
            sys.exit(-signum)

    def _summarize_pipeline(self, nodegraph, verbose=True):
        states = nodegraph.get_state_counts()

        if verbose:
            rows = [
                ("Number of nodes:", sum(states)),
                ("Number of done nodes:", states[nodegraph.DONE]),
                ("Number of runable nodes:", states[nodegraph.RUNABLE]),
                ("Number of queued nodes:", states[nodegraph.QUEUED]),
                ("Number of outdated nodes:", states[nodegraph.OUTDATED]),
                ("Number of failed nodes:", states[nodegraph.ERROR]),
            ]

            for message in padded_table(rows):
                self._logger.info(message)

        if states[nodegraph.ERROR]:
            self._logger.warning("Errors were detected while running pipeline")
        else:
            self._logger.info("Pipeline completed successfully")


def add_argument_groups(parser):
    add_scheduling_argument_group(parser)
    add_io_argument_group(parser)


def add_scheduling_argument_group(parser):
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


def add_io_argument_group(parser):
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

    return group
