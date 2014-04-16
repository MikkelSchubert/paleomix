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
from __future__ import print_function

import os
import errno
import Queue
import pickle
import signal
import logging
import collections
import multiprocessing

import pypeline.ui
import pypeline.logger

from pypeline.node import Node, MetaNode
from pypeline.nodegraph import NodeGraph, NodeGraphError
from pypeline.common.utilities import \
     safe_coerce_to_tuple, \
     fast_pickle_test
from pypeline.common.versions import \
    VersionRequirementError




class Pypeline:
    def __init__(self, config):
        self._nodes  = []
        self._config = config
        self._logger = logging.getLogger(__name__)
        # Set if a keyboard-interrupt (SIGINT) has been caught
        self._interrupted = False


    def add_nodes(self, *nodes):
        for subnodes in safe_coerce_to_tuple(nodes):
            for node in safe_coerce_to_tuple(subnodes):
                if not isinstance(node, Node):
                    raise TypeError("Node object expected, recieved %s" % repr(node))
                self._nodes.append(node)


    def run(self, max_running = 1, dry_run = False, progress_ui = "verbose"):
        try:
            nodegraph = NodeGraph(self._nodes)
        except NodeGraphError, error:
            self._logger.error(error)
            return False

        for node in nodegraph.iterflat():
            if (node.threads > max_running) and not isinstance(node, MetaNode):
                message = "Node(s) use more threads than the max allowed; " \
                          "the pipeline may therefore use more than the " \
                          "expected number of threads."
                self._logger.warning(message)
                break

        if dry_run:
            ui = pypeline.ui.VerboseUI()
            nodegraph.add_state_observer(ui)
            ui.flush()
            self._logger.info("Dry run done ...")
            return True

        old_handler = signal.signal(signal.SIGINT, self._sigint_handler)
        try:
            return self._run(nodegraph, max_running, progress_ui)
        finally:
            signal.signal(signal.SIGINT, old_handler)

        return False


    def _run(self, nodegraph, max_running, progress_ui):
        # Dictionary of nodes -> async-results
        running   = {}
        # Set of remaining nodes to be run
        remaining = set(nodegraph.iterflat())
        queue     = multiprocessing.Queue()
        pool      = multiprocessing.Pool(max_running, _init_worker, (queue,))

        errors_occured  = False
        # Signifies whether or not node-states have been refreshed following
        # the completion of the last node, in order to catch changes to the
        # states of nodes that were previously marked as DONE
        has_refreshed   = False
        # Set to true if any new nodes were started in a cycle
        has_started_any = False

        ui = pypeline.ui.get_ui(progress_ui)
        nodegraph.add_state_observer(ui)
        while running or (remaining and not self._interrupted):
            errors_occured |= not self._poll_running_nodes(running, nodegraph, queue)
            if not self._interrupted: # Prevent starting of new nodes
                if self._start_new_tasks(remaining, running, nodegraph, max_running, pool):
                    has_started_any = True
                    has_refreshed = False
                elif has_started_any and not has_refreshed:
                    # Double-check that everything is in order
                    remaining = set(nodegraph.iterflat())
                    nodegraph.refresh_states()
                    has_refreshed = True

            if running:
                ui.flush()

        pool.close()
        pool.join()

        ui.flush()
        ui.finalize()

        return not errors_occured


    def _start_new_tasks(self, remaining, running, nodegraph, max_threads, pool):
        started_nodes  = []
        idle_processes = max_threads - sum(node.threads for (node, _) in running.itervalues())
        for node in remaining:
            if not running or (idle_processes >= node.threads):
                state = nodegraph.get_node_state(node)
                if (state == nodegraph.RUNABLE):
                    try:
                        # Ensure that the node can be used in a multiprocessing context
                        fast_pickle_test(node)
                    except pickle.PicklingError, error:
                        self._logger.error("Node cannot be pickled, please file a bug-report:\n"
                                           "\tNode: %s\n\tError: %s" % (self, error))
                        nodegraph.set_node_state(node, nodegraph.ERROR)
                        started_nodes.append(node)
                        continue

                    key          = id(node)
                    proc_args    = (key, node, self._config)
                    running[key] = (node, pool.apply_async(_call_run, args = proc_args))
                    started_nodes.append(node)

                    nodegraph.set_node_state(node, nodegraph.RUNNING)
                    idle_processes -= node.threads
                elif (state in (nodegraph.DONE, nodegraph.ERROR)):
                    started_nodes.append(node)
            elif not idle_processes:
                break

        for node in started_nodes:
            remaining.remove(node)

        return bool(remaining)


    def _poll_running_nodes(self, running, nodegraph, queue):
        errors, blocking = None, True
        while running:
            node, proc = self._get_finished_node(queue, running, blocking)
            if not node:
                break

            try:
                # Re-raise exceptions from the node-process
                proc.get()
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception, errors:
                nodegraph.set_node_state(node, nodegraph.ERROR)
                self._logger.error("%s: Error occurred running command:\n%s\n",
                                   node, "\n".join(("\t" + line) for line in str(errors).strip().split("\n")))
                continue
            nodegraph.set_node_state(node, nodegraph.DONE)
            blocking = False # only block first cycle

        return not errors

    @property
    def nodes(self):
        return set(self._nodes)

    def walk_nodes(self, func):
        skip_nodes = set()

        def _walk_nodes(nodes):
            for node in nodes:
                if node in skip_nodes:
                    continue
                elif not func(node):
                    return False

                skip_nodes.add(node)
                if not _walk_nodes(node.subnodes):
                    return False
                elif not _walk_nodes(node.dependencies):
                    return False
            return True

        _walk_nodes(self._nodes)

    def list_output_files(self):
        output_files = set()

        def collect_output_files(node):
            output_files.update(node.output_files)
            return True

        self.walk_nodes(collect_output_files)

        return frozenset(map(os.path.abspath, output_files))

    def list_required_executables(self):
        requirements = collections.defaultdict(set)

        def collect_requirements(node):
            for executable in node.executables:
                _ = requirements[executable]

            for requirement in node.requirements:
                requirements[requirement.name].add(requirement)
            return True

        self.walk_nodes(collect_requirements)
        return requirements

    def print_output_files(self, print_func=print):
        for filename in sorted(self.list_output_files()):
            print_func(filename)

    def print_required_executables(self, print_func=print):
        rows = []
        pipeline_executables = self.list_required_executables()
        print_func("{: <40s} {: <10s} {}".format("Executable",
                                                 "Version",
                                                 "Required version"))

        for (name, requirements) in sorted(pipeline_executables.items()):
            if not requirements:
                print_func(name)

            for requirement in requirements:
                try:
                    version = ".".join(map(str, requirement.version))
                except VersionRequirementError:
                    version = "UNKNOWN"

                required = requirement._reqs.description
                print_func("{: <40s} {: <10s} {}".format(name,
                                                         version,
                                                         required))

    def _sigint_handler(self, signum, frame):
        """Signal handler; see signal.signal."""
        if not self._interrupted:
            self._interrupted = True
            self._logger.error("\nKeyboard interrupt detected, waiting for current tasks to complete ...\n"
                               "\t- Press CTRL-C again to force termination.\n")
        else:
            raise signal.default_int_handler(signum, frame)

    @classmethod
    def _get_finished_node(self, queue, running, blocking = True):
        """Returns a tuple containing a node that has finished running
        and it's async-result, or None for both if no such node could
        be found (and blocking is False), or if an interrupt occured
        while waiting for a node to finish.

        If blocking is True, the function will only return once a
        result becomes available, or if an interrupt occurs."""
        try:
            key = queue.get(blocking)
            return running.pop(key)
        except IOError, error:
            # User pressed ctrl-c (SIGINT), or similar event ...
            if error.errno != errno.EINTR:
                raise
        except Queue.Empty:
            pass
        return None, None



def _init_worker(queue):
    """Init function for subprocesses created by multiprocessing.Pool: Ensures that KeyboardInterrupts
    only occur in the main process, allowing us to do proper cleanup."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    # This is a workaround needed to avoid having to use multiprocessing.Manager
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
    finally:
        # See comment in _init_worker
        _call_run.queue.put(key)
