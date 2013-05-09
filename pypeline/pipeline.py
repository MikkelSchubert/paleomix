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
import pypeline.nodegraph as nodegraph
from pypeline.node import Node
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


    def run(self, max_running = 1, dry_run = False, terminate_on_error = False, collapse = True, verbose = True):
        try:
            nodes = nodegraph.NodeGraph(self._nodes)
        except nodegraph.NodeGraphError, error:
            ui.print_err(error, file = sys.stderr)
            return False

        for node in nodes.iterflat():
            if node.threads > max_running:
                ui.print_err("Node requires more threads than the maximum allowed:\n\t%s" % str(node), file = sys.stderr)
                return False


        if dry_run:
            ui.print_node_tree(nodes, collapse)
            ui.print_info("Dry run done ...", file = sys.stderr)
            return True

        try:
            running = {}
            errors = has_refreshed = has_started_any = False
            pool = multiprocessing.Pool(max_running, _init_worker)

            while not (errors and terminate_on_error):
                if not self._poll_running_nodes(running, nodes):
                    errors = True
                elif not self._start_new_tasks(running, nodes, max_running, pool):
                    if not has_refreshed and has_started_any:
                        # Double-check that everything is in order
                        nodes.refresh_states()
                        has_refreshed = True
                        continue
                    else:
                        break
                else:
                    has_started_any = True
                    has_refreshed = False

                ui.print_node_tree(nodes, collapse, verbose)

            ui.print_node_tree(nodes, collapse)
            pool.close()
            pool.join()
            
            if errors or not self._poll_running_nodes(running, nodes):
                ui.print_err("Errors were detected ...", file = sys.stderr)
                return False
        except nodegraph.NodeGraphError, errors:
            errors = "\n".join(("\t" + line) for line in str(errors).strip().split("\n"))
            ui.print_err("Error in task-graph, terminating gracefully:\n%s" % errors, file = sys.stderr)
            pool.terminate()
            pool.join()

        except KeyboardInterrupt, errors:
            ui.print_err("Keyboard interrupt detected, terminating ...", file = sys.stderr)
            pool.terminate()
            pool.join()

        ui.print_msg("Done ...", file = sys.stderr)
        return not errors


    def _start_new_tasks(self, running, nodes, max_threads, pool):
        any_runable_left = False
        idle_processes = max_threads - sum(node.threads for node in running)
        for node in nodes.iterflat():
            any_runable_left |= (nodes.get_node_state(node) in (nodes.RUNABLE, nodes.RUNNING))
            
            if idle_processes and (nodes.get_node_state(node) == nodes.RUNABLE):
                if idle_processes >= node.threads:
                    running[node] = pool.apply_async(_call_run, args = (node, self._config))
                    nodes.set_node_state(node, nodes.RUNNING)
                    idle_processes -= node.threads
            
            if any_runable_left and not idle_processes:
                break

        return any_runable_left


    @classmethod
    def _poll_running_nodes(cls, running, nodes):
        sleep_time = 0.05
        changes = errors = False
        while running and not (errors or changes):
            time.sleep(sleep_time)
            sleep_time = min(1, sleep_time * 2)

            for (node, proc) in running.items():
                if not proc.ready():
                    continue
                changes = True

                running.pop(node)                   
                try:
                    proc.get()
                except Exception, errors:
                    nodes.set_node_state(node, nodes.ERROR)                    
                    ui.print_err("%s: Error occurred running command (terminating gracefully):\n%s\n" \
                                     % (node, "\n".join(("\t" + line) for line in str(errors).strip().split("\n"))),
                                 file = sys.stderr)
                    continue
                nodes.set_node_state(node, nodes.DONE)

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
