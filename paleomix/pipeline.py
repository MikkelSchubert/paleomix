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
import collections
import errno
import logging
import multiprocessing
import os
from shlex import quote
import signal
import sys
import traceback

from queue import Empty

import paleomix.common.logging
import paleomix.common.system

from paleomix.node import Node, NodeError
from paleomix.nodegraph import FileStatusCache, NodeGraph, NodeGraphError
from paleomix.common.text import padded_table
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.versions import VersionRequirementError


class Pypeline:
    def __init__(
        self,
        nodes,
        temp_root="/tmp",
        max_threads=1,
        implicit_dependencies=False,
    ):
        if max_threads < 1:
            raise ValueError("Max threads must be >= 1")

        self._nodes = safe_coerce_to_tuple(nodes)
        for node in self._nodes:
            if not isinstance(node, Node):
                raise TypeError("Node object expected, recieved %s" % repr(node))

        self._logger = logging.getLogger(__name__)
        # Set if a keyboard-interrupt (SIGINT) has been caught
        self._interrupted = False
        # Dict of id(Node): (Node, Process)
        self._running = {}
        self._max_threads = max_threads
        self._temp_root = temp_root
        self._implicit_dependencies = implicit_dependencies

    def run(self, mode="run"):
        if mode not in ("run", "dry_run", "input_files", "output_files", "executables"):
            raise ValueError("Unknown pipeline mode {!r}".format(mode))
        elif mode == "dry_run":
            self._logger.info("Dry running pipeline ..")
        elif mode == "run":
            self._logger.info("Starting pipeline ..")
        else:
            return self._print_files(mode)

        try:
            nodegraph = NodeGraph(
                nodes=self._nodes,
                implicit_dependencies=self._implicit_dependencies,
            )
        except NodeGraphError as error:
            self._logger.error(error)
            return 1

        for node in nodegraph.iterflat():
            if node.threads > self._max_threads:
                self._logger.warn(
                    "One or more tasks require more threads than the user-defined "
                    "maximum; the number of threads used will therefore exceed the "
                    "user-defined maximum when running those tasks."
                )
                break

        if mode == "dry_run":
            self._summarize_pipeline(nodegraph)
            self._logger.info("Dry run done")
            return 0

        old_handler = signal.signal(signal.SIGINT, self._sigint_handler)

        try:
            return self._run(nodegraph)
        finally:
            signal.signal(signal.SIGINT, old_handler)
            for filename in paleomix.common.logging.get_logfiles():
                self._logger.info("Log-file written to %r", filename)

    def _run(self, nodegraph):
        # Set of remaining nodes to be run
        remaining = set(nodegraph.iterflat())

        any_errors = False
        any_changes = True
        queue = multiprocessing.Queue()
        while self._running or (remaining and not self._interrupted):
            if any_changes and (remaining and not self._interrupted):
                self._start_new_tasks(remaining, nodegraph, queue)
                if remaining and not self._running:
                    self._logger.critical("BUG: Failed to start new tasks!")
                    any_errors = True
                    break

            any_changes, any_node_errors = self._poll_running_nodes(nodegraph, queue)
            any_errors |= any_node_errors

        self._summarize_pipeline(nodegraph, verbose=any_errors)

        return 1 if any_errors else 0

    def _start_new_tasks(self, remaining, nodegraph, queue):
        idle_processes = self._max_threads - sum(
            node.threads for (node, _) in self._running.values()
        )

        # Individual node.threads may be > than max_threads, resulting in idle <= 0
        if idle_processes <= 0:
            return False

        for node in sorted(remaining, key=lambda node: node.id):
            # Any node is accepted if no nodes are running. This ensures that nodes with
            # more threads than max_threads will still be run.
            if not self._running or (idle_processes >= node.threads):
                state = nodegraph.get_node_state(node)
                if state == nodegraph.RUNABLE:
                    key = id(node)
                    proc = multiprocessing.Process(
                        target=_node_wrapper,
                        args=(queue, key, node, self._temp_root),
                        daemon=True,
                    )

                    self._running[key] = (node, proc)
                    remaining.remove(node)
                    proc.start()

                    nodegraph.set_node_state(node, nodegraph.RUNNING)
                    idle_processes -= node.threads
                elif state in (nodegraph.DONE, nodegraph.ERROR):
                    remaining.remove(node)
            elif idle_processes <= 0:
                break

    def _poll_running_nodes(self, nodegraph, queue):
        # Check for processes killed by external means to avoid waiting endlessly
        for key, (_node, proc) in self._running.items():
            proc.join(0)
            if proc.exitcode:
                message = "Process exited unexpectedly with exit code {}".format(
                    proc.exitcode
                )
                queue.put((key, NodeError(message), None))

        blocking = False
        any_errors = False
        any_changes = False
        while self._running:
            try:
                key, error, backtrace = queue.get(blocking, 1.0)
            except IOError as error:
                # User pressed ctrl-c (SIGINT), or similar event
                if error.errno != errno.EINTR:
                    raise

                continue
            except Empty:
                # Stop looping if we've already waited once or if we've finished one or
                # more nodes, since the latter means that we can maybe start new nodes
                if blocking or any_changes:
                    break

                blocking = True
                continue

            blocking = False
            any_changes = True
            node, proc = self._running.pop(key)
            proc.join()

            if error is None:
                nodegraph.set_node_state(node, nodegraph.DONE)
            else:
                nodegraph.set_node_state(node, nodegraph.ERROR)
                any_errors = True

                if not isinstance(error, NodeError):
                    error = "Unhandled exception while running {}:".format(node)

                message = str(error).split("\n")

                if backtrace:
                    message.extend("".join(backtrace).rstrip().split("\n"))

                if isinstance(error, NodeError) and error.path:
                    message.append("For more information about this error, see")
                    message.append(
                        "  " + quote(os.path.join(error.path, "pipe.errors"))
                    )

                self._logger.error("\n".join(message))

        return any_changes, any_errors

    def _print_files(self, mode):
        try:
            if mode == "input_files":
                return self._print_input_files()
            elif mode == "output_files":
                return self._print_output_files()
            elif mode == "executables":
                return self._print_required_executables()
            else:
                raise ValueError(mode)
        except BrokenPipeError:
            return 0
        except NodeGraphError as error:
            self._logger.error(error)
            return 1

    def _print_input_files(self):
        self._logger.info("Collecting and printing external input files ..")
        graph = NodeGraph(
            nodes=self._nodes,
            software_checks=False,
            implicit_dependencies=self._implicit_dependencies,
        )

        input_files = set()
        output_files = set()
        for node in graph.iterflat():
            for filename in node.input_files:
                input_files.add(os.path.abspath(filename))

            for filename in node.output_files:
                output_files.add(os.path.abspath(filename))

        for filename in sorted(input_files - output_files):
            print(filename)

        return 0

    def _print_output_files(self, file=sys.stdout):
        self._logger.info("Collecting and printing output files ..")
        output_files = {}
        cache = FileStatusCache()
        graph = NodeGraph(
            nodes=self._nodes,
            software_checks=False,
            implicit_dependencies=self._implicit_dependencies,
            cache_factory=lambda: cache,
        )

        def _set_output_file_state(filenames, state):
            for filename in filenames:
                output_files[os.path.abspath(filename)] = state

        for node in graph.iterflat():
            state = graph.get_node_state(node)
            if state == NodeGraph.DONE:
                _set_output_file_state(node.output_files, "Ready      ")
                continue

            # Pending/queued nodes may have outdated output files
            missing_files = frozenset(cache.missing_files(node.output_files))
            _set_output_file_state(missing_files, "Missing    ")
            _set_output_file_state(node.output_files - missing_files, "Outdated   ")

        for filename, state in sorted(output_files.items()):
            print(state, filename, file=file)

        return 0

    def _print_required_executables(self, file=sys.stdout):
        self._logger.info("Collecting and printing required executables ..")
        graph = NodeGraph(
            nodes=self._nodes,
            software_checks=False,
            implicit_dependencies=self._implicit_dependencies,
        )

        executables = set()
        requirements = collections.defaultdict(set)
        requirement_executables = set()
        for node in graph.iterflat():
            executables.update(node.executables)

            for requirement in node.requirements:
                requirements[requirement.name].add(requirement)
                requirement_executables.add(requirement.executable)

        for executable in executables - requirement_executables:
            requirements[executable] = set()

        template = "{: <40s} {: <11s} {}"
        print(template.format("Executable", "Version", "Required version"), file=file)

        for (name, requirements) in sorted(requirements.items()):
            if not requirements:
                print(template.format(name, "-", "any version"), file=file)
                continue

            for requirement in requirements:
                try:
                    version = requirement.version_str
                except VersionRequirementError:
                    version = "ERROR"

                print(template.format(name, version, requirement.checks), file=file)

        return 0

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
            for _node, proc in self._running.values():
                proc.terminate()
            for _node, proc in self._running.values():
                proc.join()
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


def _node_wrapper(queue, key, node, temp_root):
    name = "paleomix task"
    if len(sys.argv) > 1:
        name = "paleomix {} task".format(sys.argv[1])

    paleomix.common.system.set_procname(name)
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    try:
        node.run(temp_root)
        queue.put((key, None, None))
    except NodeError as error:
        backtrace = []
        if error.__cause__ and error.__cause__.__traceback__:
            backtrace = traceback.format_tb(error.__cause__.__traceback__)
            backtrace.append("  {!r}".format(error.__cause__))

        queue.put((key, error, backtrace))
    except Exception:
        _exc_type, exc_value, exc_traceback = sys.exc_info()
        backtrace = traceback.format_tb(exc_traceback)
        backtrace.append("  {!r}".format(exc_value))

        queue.put((key, exc_value, backtrace))
