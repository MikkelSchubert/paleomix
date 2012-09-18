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


    def run(self, max_running = 1, dry_run = False, terminate_on_error = False):
        try:
            nodes = nodegraph.NodeGraph(self._nodes)
        except nodegraph.NodeGraphError, error:
            ui.print_err(error)
            return False

        for node in nodes.iterflat():
            if node.threads > max_running:
                ui.print_err("Node requires more threads than the maximum allowed:\n\t%s" % str(node))
                return 1


        if dry_run:
            ui.print_node_tree(nodes)
            ui.print_msg("Dry run done ...")
            return 0
    
        try:
            running = {}
            pool = multiprocessing.Pool(max_running, _init_worker)
            while self._poll_running_nodes(running, nodes) or not terminate_on_error:
                if not self._start_new_tasks(running, nodes, max_running, pool):
                    ui.print_node_tree(nodes)
                    break

                ui.print_node_tree(nodes)

            pool.close()
            pool.join()

            if not self._poll_running_nodes(running, nodes):
                ui.print_err("Errors were detected ...")
                return False

        except nodegraph.NodeGraphError, errors:
            errors = "\n".join(("\t" + line) for line in str(errors).strip().split("\n"))
            ui.print_err("Error in task-graph, terminating gracefully:\n%s" % errors)
            pool.terminate()
            pool.join()

        except KeyboardInterrupt:
            ui.print_err("Keyboard interrupt detected, terminating ...")
            pool.terminate()
            pool.join()

        ui.print_msg("Done ...")
        return True


    def _start_new_tasks(self, running, nodes, max_threads, pool):
        any_runable_left = False
        idle_processes = max_threads - sum(node.threads for node in running)
        for node in nodes.iterflat():
            any_runable_left |= (node.state in (node.RUNABLE, node.RUNNING))
            
            if idle_processes and (node.state == node.RUNABLE):
                if idle_processes >= node.threads:
                    running[node] = pool.apply_async(_call_run, args = (node.task, self._config))
                    nodes.set_task_state(node, node.RUNNING)
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
            sleep_time = min(5, sleep_time * 2)

            for (node, proc) in running.items():
                if not proc.ready():
                    continue
                changes = True

                running.pop(node)
                nodes.set_task_state(node, None)
                   
                try:
                    proc.get()
                except Exception, errors:
                    nodes.set_task_state(node, node.ERROR)                    
                    ui.print_err("%s: Error occurred running command (terminating gracefully):\n%s\n" \
                                     % (node, "\n".join(("\t" + line) for line in str(errors).strip().split("\n"))))

        return not errors
 



def _init_worker():
    """Init function for subprocesses created by multiprocessing.Pool: Ensures that KeyboardInterrupts 
    only occur in the main process, allowing us to do proper cleanup."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _call_run(node, config):
    """Wrapper function, required in order to call Node.run()
    in subprocesses, since it is not possible to pickle 
    bound functions (e.g. self.run)"""
    return node.run(config)
