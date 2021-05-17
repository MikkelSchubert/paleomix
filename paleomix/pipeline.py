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
import signal
import sys
import traceback

from queue import Empty

import paleomix.common.logging

from paleomix.node import Node, NodeError, NodeUnhandledException
from paleomix.nodegraph import FileStatusCache, NodeGraph, NodeGraphError
from paleomix.common.text import padded_table
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.versions import VersionRequirementError


class Pypeline:
    def __init__(self, config, implicit_dependencies=False):
        self._nodes = []
        self._config = config
        self._logger = logging.getLogger(__name__)
        # Set if a keyboard-interrupt (SIGINT) has been caught
        self._interrupted = False
        self._queue = multiprocessing.Queue()
        self._pool = None
        self._implicit_dependencies = implicit_dependencies

    def add_nodes(self, *nodes):
        for subnodes in safe_coerce_to_tuple(nodes):
            for node in safe_coerce_to_tuple(subnodes):
                if not isinstance(node, Node):
                    raise TypeError("Node object expected, recieved %s" % repr(node))
                self._nodes.append(node)

    def run(self, max_threads=1, mode="run"):
        if mode not in ("run", "dry_run", "input_files", "output_files", "executables"):
            raise ValueError("Unknown pipeline mode {!r}".format(mode))
        elif max_threads < 1:
            raise ValueError("Max threads must be >= 1")
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
            if node.threads > max_threads:
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

        self._pool = multiprocessing.Pool(max_threads, _init_worker, (self._queue,))
        old_handler = signal.signal(signal.SIGINT, self._sigint_handler)

        try:
            return self._run(nodegraph, max_threads)
        finally:
            signal.signal(signal.SIGINT, old_handler)
            for filename in paleomix.common.logging.get_logfiles():
                self._logger.info("Log-file written to %r", filename)

    def _run(self, nodegraph, max_threads):
        # Dictionary of nodes -> async-results
        running = {}
        # Set of remaining nodes to be run
        remaining = set(nodegraph.iterflat())

        is_ok = True
        while running or (remaining and not self._interrupted):
            is_ok &= self._poll_running_nodes(running, nodegraph, self._queue)

            if not self._interrupted:  # Prevent starting of new nodes
                self._start_new_tasks(
                    remaining, running, nodegraph, max_threads, self._pool
                )

        self._pool.close()
        self._pool.join()

        self._summarize_pipeline(nodegraph, verbose=not is_ok)

        return 0 if is_ok else 1

    def _start_new_tasks(self, remaining, running, nodegraph, max_threads, pool):
        started_nodes = []
        idle_processes = max_threads - sum(
            node.threads for (node, _) in running.values()
        )

        if not idle_processes:
            return False

        for node in sorted(remaining, key=lambda node: node.id):
            if not running or (idle_processes >= node.threads):
                state = nodegraph.get_node_state(node)
                if state == nodegraph.RUNABLE:
                    key = id(node)
                    proc_args = (key, node, self._config)
                    running[key] = (node, pool.apply_async(_call_run, args=proc_args))
                    started_nodes.append(node)

                    nodegraph.set_node_state(node, nodegraph.RUNNING)
                    idle_processes -= node.threads
                elif state in (nodegraph.DONE, nodegraph.ERROR):
                    started_nodes.append(node)
            elif idle_processes <= 0:
                break

        for node in started_nodes:
            remaining.remove(node)

    def _poll_running_nodes(self, running, nodegraph, queue):
        error_happened = False
        blocking = False

        while running and not error_happened:
            node, proc = self._get_finished_node(queue, running, blocking)
            if not node:
                if blocking:
                    break

                blocking = True
                continue

            try:
                # Re-raise exceptions from the node-process
                proc.get()
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as errors:
                error_happened = True
                nodegraph.set_node_state(node, nodegraph.ERROR)

                self._logger.error("%s while %s:", type(errors).__name__, node)
                for line in str(errors).strip().split("\n"):
                    self._logger.error("    %s", line)

            if not error_happened:
                nodegraph.set_node_state(node, nodegraph.DONE)

        return not error_happened

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
        graph = NodeGraph(nodes=self._nodes, software_checks=False)

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
        cache = FileStatusCache()
        graph = NodeGraph(
            nodes=self._nodes,
            software_checks=False,
            cache_factory=lambda: cache,
        )

        output_files = {}
        for node in graph.iterflat():
            state = graph.get_node_state(node)
            if state == NodeGraph.DONE:
                state = "Ready      "
            elif state == graph.RUNABLE:
                # Pending nodes may have outdated output files
                if cache.missing_files(node.output_files):
                    state = "Missing    "
                else:
                    state = "Outdated   "
            else:
                state = "Missing    "

            for filename in node.output_files:
                output_files[os.path.abspath(filename)] = state

        for filename, state in sorted(output_files.items()):
            print(state, filename, file=file)

        return 0

    def _print_required_executables(self, file=sys.stdout):
        self._logger.info("Collecting and printing required executables ..")
        graph = NodeGraph(
            nodes=self._nodes,
            software_checks=False,
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
                    if requirement.version:
                        version = "v" + ".".join(map(str, requirement.version))
                    else:
                        version = "NA"
                except VersionRequirementError:
                    version = "UNKNOWN"

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
            self._pool.terminate()
            sys.exit(1)

    @classmethod
    def _get_finished_node(cls, queue, running, blocking):
        """Returns a tuple containing a node that has finished running
        and it's async-result, or None for both if no such node could
        be found (and blocking is False), or if an interrupt occured
        while waiting for a node to finish.

        If blocking is True, the function will timeout after 0.1s.
        """
        try:
            key = queue.get(blocking, 0.1)
            return running.pop(key)
        except IOError as error:
            # User pressed ctrl-c (SIGINT), or similar event
            if error.errno != errno.EINTR:
                raise
        except Empty:
            pass
        return None, None

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


def _init_worker(queue):
    """Init function for subprocesses created by multiprocessing.Pool: Ensures
    that KeyboardInterrupts only occur in the main process, allowing us to do
    proper cleanup.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    # This is a workaround to avoid having to use multiprocessing.Manager
    # to create the Queue objects; this is needed because the Manager class
    # creates it own process, which inherits the signal-handlers of the main
    # process, causing some rather odd behavior when the user causes a SIGINT.
    _call_run.queue = queue


def _call_run(key, node, config):
    """Wrapper function, required in order to call Node.run()
    in subprocesses, since it is not possible to pickle
    bound functions (e.g. self.run)"""
    try:
        return node.run(config)
    except NodeError:
        raise
    except Exception:
        message = "Unhandled error running Node:\n\n%s" % (traceback.format_exc(),)

        raise NodeUnhandledException(message)
    finally:
        # See comment in _init_worker
        _call_run.queue.put(key)
