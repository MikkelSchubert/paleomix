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
from __future__ import print_function

import errno
import logging
import multiprocessing
import os
import pickle
import Queue
import signal
import traceback

import paleomix.ui
import paleomix.logger

from paleomix.node import \
    Node, \
    NodeError, \
    NodeUnhandledException
from paleomix.nodegraph import \
    FileStatusCache, \
    NodeGraph, \
    NodeGraphError
from paleomix.common.utilities import \
    safe_coerce_to_tuple, \
    fast_pickle_test
from paleomix.common.versions import \
    VersionRequirementError


class Pypeline(object):
    def __init__(self, config):
        self._nodes = []
        self._config = config
        self._logger = logging.getLogger(__name__)
        # Set if a keyboard-interrupt (SIGINT) has been caught
        self._interrupted = False
        self._queue = multiprocessing.Queue()
        self._pool = multiprocessing.Pool(1, _init_worker, (self._queue,))

    def add_nodes(self, *nodes):
        for subnodes in safe_coerce_to_tuple(nodes):
            for node in safe_coerce_to_tuple(subnodes):
                if not isinstance(node, Node):
                    raise TypeError("Node object expected, recieved %s"
                                    % repr(node))
                self._nodes.append(node)

    def run(self, max_threads=1, dry_run=False, progress_ui="verbose"):
        if max_threads < 1:
            raise ValueError("Max threads must be >= 1")
        _update_nprocesses(self._pool, max_threads)

        try:
            nodegraph = NodeGraph(self._nodes)
        except NodeGraphError, error:
            self._logger.error(error)
            return False

        for node in nodegraph.iterflat():
            if (node.threads > max_threads):
                message = "Node(s) use more threads than the max allowed; " \
                          "the pipeline may therefore use more than the " \
                          "expected number of threads.\n"
                paleomix.ui.print_warn(message)
                break

        if dry_run:
            progress_printer = paleomix.ui.QuietUI()
            nodegraph.add_state_observer(progress_printer)
            progress_printer.flush()
            progress_printer.finalize()
            self._logger.info("Dry run done ...")
            return True

        old_handler = signal.signal(signal.SIGINT, self._sigint_handler)
        try:
            return self._run(nodegraph, max_threads, progress_ui)
        finally:
            signal.signal(signal.SIGINT, old_handler)

        return False

    def _run(self, nodegraph, max_threads, progress_ui):
        # Dictionary of nodes -> async-results
        running = {}
        # Set of remaining nodes to be run
        remaining = set(nodegraph.iterflat())

        is_ok = True
        progress_printer = paleomix.ui.get_ui(progress_ui)
        progress_printer.max_threads = max_threads
        nodegraph.add_state_observer(progress_printer)

        with paleomix.ui.CommandLine() as cli:
            while running or (remaining and not self._interrupted):
                is_ok &= self._poll_running_nodes(running,
                                                  nodegraph,
                                                  self._queue)

                if not self._interrupted:  # Prevent starting of new nodes
                    self._start_new_tasks(remaining, running, nodegraph,
                                          max_threads, self._pool)

                if running:
                    progress_printer.flush()

                max_threads = cli.process_key_presses(nodegraph,
                                                      max_threads,
                                                      progress_printer)
                progress_printer.max_threads = max_threads
                _update_nprocesses(self._pool, max_threads)

        self._pool.close()
        self._pool.join()

        progress_printer.flush()
        progress_printer.finalize()

        return is_ok

    def _start_new_tasks(self, remaining, running, nodegraph, max_threads,
                         pool):
        started_nodes = []
        idle_processes = max_threads \
            - sum(node.threads for (node, _) in running.itervalues())

        if not idle_processes:
            return False

        for node in remaining:
            if not running or (idle_processes >= node.threads):
                state = nodegraph.get_node_state(node)
                if state == nodegraph.RUNABLE:
                    try:
                        # The multi-processing module relies on pickling
                        fast_pickle_test(node)
                    except pickle.PicklingError, error:
                        self._logger.error("Node cannot be pickled; please "
                                           "file a bug-report:\n"
                                           "\tNode: %s\n\tError: %s"
                                           % (self, error))
                        nodegraph.set_node_state(node, nodegraph.ERROR)
                        started_nodes.append(node)
                        continue

                    key = id(node)
                    proc_args = (key, node, self._config)
                    running[key] = (node, pool.apply_async(_call_run,
                                                           args=proc_args))
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
        errors = None
        blocking = False

        while running and not errors:
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
            except Exception, errors:
                nodegraph.set_node_state(node, nodegraph.ERROR)

                message = [str(node),
                           "  Error (%r) occurred running command:"
                           % (type(errors).__name__)]

                for line in str(errors).strip().split("\n"):
                    message.append("    %s" % (line,))
                message.append("")

                self._logger.error("\n".join(message))

            if not errors:
                nodegraph.set_node_state(node, nodegraph.DONE)

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
                if not _walk_nodes(node.dependencies):
                    return False
            return True

        _walk_nodes(self._nodes)

    def list_input_files(self):
        """Returns a set containing the absolute path of all input files
        required by the current pipeline. These do not include any file
        generated by the pipeline itself (output files).
        """
        input_files = set()
        output_files = set()

        def collect_output_files(node):
            for filename in node.input_files:
                input_files.add(os.path.abspath(filename))

            for filename in node.output_files:
                output_files.add(os.path.abspath(filename))

            return True

        self.walk_nodes(collect_output_files)

        return input_files - output_files

    def list_output_files(self):
        cache = FileStatusCache()
        nodegraph = NodeGraph(self._nodes, lambda: cache)
        output_files = {}

        def collect_output_files(node):
            state = None
            if nodegraph.is_done(node, cache):
                state = nodegraph.DONE
                if nodegraph.is_outdated(node, cache):
                    state = nodegraph.OUTDATED

            for filename in node.output_files:
                output_files[os.path.abspath(filename)] = state

            return True

        self.walk_nodes(collect_output_files)

        return output_files

    def list_required_executables(self):
        requirements = {}

        def collect_requirements(node):
            for executable in node.executables:
                if executable not in requirements:
                    requirements[executable] = set()

            for requirement in node.requirements:
                if requirement.name not in requirements:
                    requirements[requirement.name] = set()

                requirements[requirement.name].add(requirement)

                executable = requirement.executable
                if not requirements.get(executable):
                    requirements.pop(executable, None)

            return True

        self.walk_nodes(collect_requirements)
        return requirements

    def print_output_files(self, print_func=print):
        output_files = self.list_output_files()

        for filename, state in sorted(output_files.iteritems()):
            if state == NodeGraph.DONE:
                state = "Ready      "
            elif state == NodeGraph.OUTDATED:
                state = "Outdated   "
            else:
                state = "Missing    "

            print_func("%s\t%s" % (state, filename))

    def print_input_files(self, print_func=print):
        """Prints the absolute path of all input files required by the current
        pipeline, excluding any file generated by the pipeline itself (output
        files). One file is printed per line.
        """
        input_files = self.list_input_files()
        for filename in sorted(input_files):
            print_func("%s" % (filename,))

    def print_required_executables(self, print_func=print):
        template = "{: <40s} {: <11s} {}"
        pipeline_executables = self.list_required_executables()
        print_func(template.format("Executable",
                                   "Version",
                                   "Required version"))

        for (name, requirements) in sorted(pipeline_executables.items()):
            if not requirements:
                print_func(template.format(name, "-", "any version"))
                continue

            for requirement in requirements:
                try:
                    if requirement.version:
                        version = "v" + ".".join(map(str, requirement.version))
                    else:
                        version = "NA"
                except VersionRequirementError:
                    version = "UNKNOWN"

                print_func(template.format(name, version, requirement.checks))

    def _sigint_handler(self, signum, frame):
        """Signal handler; see signal.signal."""
        if not self._interrupted:
            self._interrupted = True
            self._logger.error("\nKeyboard interrupt detected, waiting for "
                               "running tasks to complete ... Press CTRL-C "
                               "again to force termination.\n")
        else:
            self._pool.terminate()
            raise signal.default_int_handler(signum, frame)

    def to_dot(self, destination):
        """Writes a simlpe dot file to the specified destination, representing
        the full dependency tree. Nodes are named by their class.
        """
        try:
            nodegraph = NodeGraph(self._nodes)
        except NodeGraphError, error:
            self._logger.error(error)
            return False

        # Dict recording all dependencies of nodes
        meta_dependencies = {}
        # Dict recording if anything depends on a speific node
        meta_rev_dependencies = {}
        for node in nodegraph.iterflat():
            selection = node.dependencies
            meta_dependencies[node] = selection
            for dep in selection:
                meta_rev_dependencies[dep] = True

        return self._write_dot(destination,
                               meta_dependencies,
                               meta_rev_dependencies)

    @classmethod
    def _write_dot(cls, destination, meta_dependencies, meta_rev_dependencies):
        """Writes simple dot file, in which each node is connected to their
        dependencies, using the object IDs as the node names. Labels are
        derived from the class names, excluding any "Node" postfix.
        """
        with open(destination, "w") as out:
            out.write("digraph G {\n")
            out.write("  graph [ dpi = 75 ];\n")
            out.write("  node [shape=record,width=.1,height=.1];\n")
            out.write("  splines=ortho;\n\n")

            for node, dependencies in meta_dependencies.iteritems():
                node_id = "Node_%i" % (id(node),)
                node_type = node.__class__.__name__
                if node_type.endswith("Node"):
                    node_type = node_type[:-4]

                rank = None
                color = "white"
                if not meta_dependencies.get(node):
                    color = "red"
                elif not meta_rev_dependencies.get(node):
                    color = "green"
                    rank = "sink"

                if rank is not None:
                    out.write("  {")
                    out.write("    rank = %s;\n  " % (rank,))

                out.write('  %s [label="%s"; fillcolor=%s; style=filled]\n'
                          % (node_id, node_type, color))

                if rank is not None:
                    out.write("  }")

                for dependency in dependencies:
                    dep_id = "Node_%i" % (id(dependency),)
                    out.write("  %s -> %s\n" % (dep_id, node_id))
                out.write("\n")

            out.write("}\n")

        return True

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
        except IOError, error:
            # User pressed ctrl-c (SIGINT), or similar event ...
            if error.errno != errno.EINTR:
                raise
        except Queue.Empty:
            pass
        return None, None


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
        message = "Unhandled error running Node:\n\n%s" \
            % (traceback.format_exc(),)

        raise NodeUnhandledException(message)
    finally:
        # See comment in _init_worker
        _call_run.queue.put(key)


def _update_nprocesses(pool, processes):
    """multiprocessing.Pool does not expose calls to change number of active
    processes, but does in fact support this for the 'maxtasksperchild' option.
    This function calls the related private functions to increase the number
    of available processes."""
    # FIXME: Catch ERRNO 11:
    # OSError: [Errno 11] Resource temporarily unavailable
    if pool._processes < processes:
        pool._processes = processes
        pool._repopulate_pool()
