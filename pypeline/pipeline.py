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
import sys
import time
import signal
import multiprocessing

import pypeline.ui as ui
from pypeline.node import Node
from pypeline.nodegraph import NodeGraph, NodeGraphError
from pypeline.common.utilities import safe_coerce_to_tuple




class Pypeline:
    def __init__(self, config):
        self._nodes = []
        self._config = config


    def add_nodes(self, *nodes):
        for subnodes in safe_coerce_to_tuple(nodes):
            for node in safe_coerce_to_tuple(subnodes):
                if not isinstance(node, Node):
                    raise TypeError("Node object expected, recieved %s" % repr(node))
                self._nodes.append(node)


    def run(self, max_running = 1, dry_run = False, collapse = True, verbose = True):
        try:
            nodegraph = NodeGraph(self._nodes)
        except NodeGraphError, error:
            ui.print_err(error, file = sys.stderr)
            return False

        remaining = set(nodegraph.iterflat())
        for node in remaining:
            if node.threads > max_running:
                ui.print_err("Node requires more threads than the maximum allowed:\n\t%s" \
                             % str(node), file = sys.stderr)
                return False

        if dry_run:
            ui.print_node_tree(nodegraph, collapse)
            ui.print_info("Dry run done ...", file = sys.stderr)
            return True

        running = {}
        interrupted_once = errors = has_refreshed = has_started_any = False
        pool = multiprocessing.Pool(max_running, _init_worker)

        while running or remaining:
            try:
                errors |= not self._poll_running_nodes(running, nodegraph)
                if not interrupted_once: # Prevent starting of new nodes
                    if self._start_new_tasks(remaining, running, nodegraph, max_running, pool):
                        has_started_any = True
                        has_refreshed = False
                    elif has_started_any and not has_refreshed:
                        # Double-check that everything is in order
                        remaining = set(nodegraph.iterflat())
                        nodegraph.refresh_states()
                        has_refreshed = True

                if running:
                    ui.print_node_tree(nodegraph, collapse, verbose)
            except KeyboardInterrupt:
                if interrupted_once:
                    ui.print_err("\nTerminating now!\n", file = sys.stderr)
                    pool.terminate()
                    pool.join()
                    return False

                remaining, interrupted_once = set(), True
                ui.print_err("\nKeyboard interrupt detected, waiting for current tasks to complete ...",
                             file = sys.stderr)
                ui.print_err("\t- Press CTRL-C again to force termination.\n",
                             file = sys.stderr)

        ui.print_node_tree(nodegraph, collapse)
        pool.close()
        pool.join()

        if errors:
            ui.print_err("Errors were detected ...", file = sys.stderr)
        ui.print_msg("Done ...", file = sys.stderr)
        return not errors


    def _start_new_tasks(self, remaining, running, nodegraph, max_threads, pool):
        started_nodes  = []
        idle_processes = max_threads - sum(node.threads for node in running)
        for node in remaining:
            if (idle_processes >= node.threads):
                state = nodegraph.get_node_state(node)
                if (state == nodegraph.RUNABLE):
                    started_nodes.append(node)
                    running[node] = pool.apply_async(_call_run, args = (node, self._config))
                    nodegraph.set_node_state(node, nodegraph.RUNNING)
                    idle_processes -= node.threads
                elif (state in (nodegraph.DONE, nodegraph.ERROR)):
                    started_nodes.append(node)
            elif not idle_processes:
                break

        for node in started_nodes:
            remaining.remove(node)

        return bool(remaining)


    @classmethod
    def _poll_running_nodes(cls, running, nodegraph):
        sleep_time = 0.05
        changes = errors = False
        while running and not (errors or changes):
            for (node, proc) in running.items():
                if not proc.ready():
                    continue
                changes = True

                try:
                    proc.get()
                except (KeyboardInterrupt, SystemExit):
                    raise
                except Exception, errors:
                    nodegraph.set_node_state(node, nodegraph.ERROR)
                    running.pop(node)
                    ui.print_err("%s: Error occurred running command:\n%s\n" \
                                     % (node, "\n".join(("\t" + line) for line in str(errors).strip().split("\n"))),
                                 file = sys.stderr)
                    continue
                nodegraph.set_node_state(node, nodegraph.DONE)
                running.pop(node)

            if not (errors or changes):
                time.sleep(sleep_time)
                sleep_time = min(1, sleep_time * 2)

        return not errors

    @property
    def nodes(self):
        return set(self._nodes)



def _init_worker():
    """Init function for subprocesses created by multiprocessing.Pool: Ensures that KeyboardInterrupts
    only occur in the main process, allowing us to do proper cleanup."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _call_run(node, config):
    """Wrapper function, required in order to call Node.run()
    in subprocesses, since it is not possible to pickle
    bound functions (e.g. self.run)"""
    return node.run(config)
