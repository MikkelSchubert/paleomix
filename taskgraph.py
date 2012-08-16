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
import os
import collections

from pypeline.common.fileutils import missing_executables



class TaskError(RuntimeError):
    pass
    


class TaskGraph:
    class Node:
        DONE, RUNNING, RUNABLE, QUEUED, OUTDATED, ERROR = range(6)

        def __init__(self, task):
            self.task          = task
            self._state        = None
            self._fixed_state  = None
            self.subnodes      = frozenset()
            self.dependencies  = frozenset()
            self.threads       = int(task.threads)


        @property
        def state(self):
            if self._fixed_state is not None:
                return self._fixed_state

            assert self._state is not None
            return self._state


        def __str__(self):
            return str(self.task)


    def __init__(self, tasks):
        self._tasks = {}
        def collapse(subtasks):
            for task in subtasks:
                # TODO: Check that all are 'Node's
                if task not in self._tasks:
                    self._tasks[task] = TaskGraph.Node(task)
                    collapse(task.subnodes)
                    collapse(task.dependencies)
        collapse(tasks)

        self._graph = self._build_graph(self._tasks)
        self._graph_valid = False

        self._check_file_dependencies(self._tasks)
        self._check_required_executables(self._tasks)


    def set_task_state(self, node, state):
        if state not in (None, TaskGraph.Node.RUNNING, TaskGraph.Node.ERROR):
            raise ValueError("Cannot set states other than RUNNING and ERROR, or cleared (None).")
        
        self._tasks[node.task]._fixed_state = state
        self._graph_valid = False


    def __iter__(self):
        """Returns a graph of tasks."""
        self._update_graph()
        return iter(self._graph)


    def iterflat(self):
        self._update_graph()
        return iter(self._tasks.values())


    def _update_graph(self):
        if not self._graph_valid:
            for node in self._tasks.itervalues():
                node._state = None

            for node in self._graph:
                self._update_states(node)

            self._graph_valid = True
        

    @classmethod
    def _build_graph(cls, tasks):
        cls._check_file_dependencies(tasks)

        for (task, node) in tasks.iteritems():
            node.subnodes     = frozenset(tasks[subnode] for subnode in task.subnodes)
            node.dependencies = frozenset(tasks[dependency] for dependency in task.dependencies)

        return cls._get_top_nodes(tasks)


    @classmethod
    def _get_top_nodes(cls, tasks):
        """Returns a sequence of nodes that are not depended upon by any other nodes."""

        dependencies = set()
        for (_, node) in tasks.iteritems():
            dependencies.update(node.subnodes)
            dependencies.update(node.dependencies)
        
        return set(tasks.itervalues()) - dependencies
        
 
    @classmethod
    def _update_states(cls, node):
        if node._state is not None:
            # Possibly return a fixed state
            return node.state

        # Update sub-tasks, before checking for fixed states
        state = node.DONE
        for subnode in node.subnodes:
            state = max(state, cls._update_states(subnode))
        for subnode in node.dependencies:
            state = max(state, cls._update_states(subnode))

        if state == node.DONE:
            if not node.task.is_done or node.task.is_outdated:
                state = node.RUNABLE
        elif state in (node.RUNNING, node.RUNABLE, node.QUEUED):
            if node.task.is_done:
                state = node.OUTDATED
            else:
                state = node.QUEUED
        node._state = state

        # Possibly return a fixed state
        return node.state


    @classmethod
    def _check_required_executables(cls, tasks):
        missing_exec = set()
        for task in tasks:
            missing_exec.update(missing_executables(task.executables))


        if missing_exec:
            raise TaskError("Required executables are missing:\n\t%s" \
                                % ("\n\t".join(sorted(missing_exec))))
            

    @classmethod
    def _check_file_dependencies(cls, tasks):
        input_files = collections.defaultdict(list)
        output_files = collections.defaultdict(list)

        for task in tasks:
            for filename in task.input_files:
                input_files[filename].append(task)
            
            for filename in task.output_files:
                output_files[filename].append(task)

        error_messages = []
        error_messages.extend(cls._check_output_files(output_files))
        error_messages.extend(cls._check_input_dependencies(input_files, output_files))

        if error_messages:
            messages = []
            for error in error_messages:
                for line in error.split("\n"):
                    messages.append("\t" + line)

            raise TaskError("Errors detected during graph construction:\n%s" \
                                % ("\n".join(messages)),)

    @classmethod
    def _check_output_files(cls, output_files):
        for (filename, tasks) in output_files.iteritems():
            if (len(tasks) > 1):
                yield "%i nodes clobber a file: %s:\n\t%s" \
                    % (len(tasks), filename, "\n\t".join(str(task) for task in tasks))

    @classmethod
    def _check_input_dependencies(cls, input_files, output_files):
        for (filename, nodes) in input_files.iteritems():
            if (filename in output_files):
                producer = output_files[filename][0]
                for consumer in nodes:
                    if not cls._node_has_dependency(consumer, producer):
                        yield "Node depends on dynamically created file, but not on the node creating it:" + \
                            "\n\tDependent node: %s\n\tFilename: %s\n\tCreated by: %s" \
                            % (consumer, filename, producer)
            elif not os.path.exists(filename):
                yield "Required file does not exist, and is not created by a task: %s" % (filename, )


    @classmethod
    def _node_has_dependency(cls, node, dependency):
        subnodes = node.dependencies | node.subnodes
        if dependency in subnodes:
            return True

        for subnode in subnodes:
            if cls._node_has_dependency(subnode, dependency):
                return True

        return False
