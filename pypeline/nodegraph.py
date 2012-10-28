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

import ui
from pypeline.common.fileutils import missing_executables


# Max number of error messages of each type
_MAX_ERROR_MESSAGES = 10


class NodeGraphError(RuntimeError):
    pass
    


class NodeGraph:
    DONE, RUNNING, RUNABLE, QUEUED, OUTDATED, ERROR = range(6)

    def __init__(self, nodes):
        self._reverse_dependencies = collections.defaultdict(set)
        self._collect_reverse_dependencies(nodes, self._reverse_dependencies)
        self._top_nodes = [node for (node, rev_deps) in self._reverse_dependencies.iteritems() if not rev_deps]

        self._check_file_dependencies(self._reverse_dependencies)
        self._check_required_executables(self._reverse_dependencies)

        self._graph_valid = False
        self._fixed_states = {}
        self._states = {}


    def get_node_state(self, node):
        if self._fixed_states.get(node) is not None:
            return self._fixed_states[node]
        return self._states.get(node)


    def set_node_state(self, node, state):
        if state not in (None, NodeGraph.RUNNING, NodeGraph.ERROR):
            raise ValueError("Cannot set states other than RUNNING and ERROR, or cleared (None).")
        
        self._fixed_states[node] = state
        self._graph_valid = False


    def __iter__(self):
        """Returns a graph of nodes."""
        self._update_graph()
        return iter(self._top_nodes)


    def iterflat(self):
        self._update_graph()
        return iter(self._reverse_dependencies)


    def _update_graph(self):
        if not self._graph_valid:
            self._states = {}
            for node in self._reverse_dependencies:
                self._update_states(node)

            self._graph_valid = True
        

    def _update_states(self, node):
        state = self.get_node_state(node)
        if state is not None:
            # Possibly return a fixed state
            return state
 
        # Update sub-nodes before checking for fixed states
        state = NodeGraph.DONE
        for subnode in (node.subnodes | node.dependencies):
            state = max(state, self._update_states(subnode))

        try:
            if state == NodeGraph.DONE:
                if not node.is_done or node.is_outdated:
                    state = NodeGraph.RUNABLE
            elif state in (NodeGraph.RUNNING, NodeGraph.RUNABLE, NodeGraph.QUEUED):
                if node.is_done:
                    state = NodeGraph.OUTDATED
                else:
                    state = NodeGraph.QUEUED
        except OSError, e:
            # Typically hapens if base input files are removed, causing a node that
            # 'is_done' to call modified_after on missing files in 'is_outdated'
            ui.print_err("OSError checking state of Node: %s" % e)
            state = NodeGraph.ERROR
        self._states[node] = state
 
        # Possibly return a fixed state
        return self.get_node_state(node)


    @classmethod
    def _check_required_executables(cls, nodes):
        missing_exec = set()
        for node in nodes:
            missing_exec.update(missing_executables(node.executables))

        if missing_exec:
            raise NodeGraphError("Required executables are missing:\n\t%s" \
                                % ("\n\t".join(sorted(missing_exec))))
            

    @classmethod
    def _check_file_dependencies(cls, nodes):
        input_files = collections.defaultdict(list)
        output_files = collections.defaultdict(list)

        for node in nodes:
            for filename in node.input_files:
                input_files[filename].append(node)
            
            for filename in node.output_files:
                output_files[filename].append(node)

        max_messages = range(_MAX_ERROR_MESSAGES)
        error_messages = []
        error_messages.extend(zip(max_messages, cls._check_output_files(output_files)))
        error_messages.extend(zip(max_messages, cls._check_input_dependencies(input_files, output_files, nodes)))

        if error_messages:
            messages = []
            for (_, error) in error_messages:
                for line in error.split("\n"):
                    messages.append("\t" + line)

            raise NodeGraphError("Errors detected during graph construction (max %i shown):\n%s" \
                                % (_MAX_ERROR_MESSAGES * 2, "\n".join(messages)),)

    @classmethod
    def _check_output_files(cls, output_files):
        for (filename, nodes) in output_files.iteritems():
            if (len(nodes) > 1):
                yield "%i nodes clobber a file: %s:\n\t%s" \
                    % (len(nodes), filename, "\n\t".join(str(node) for node in nodes))

    @classmethod
    def _check_input_dependencies(cls, input_files, output_files, nodes):
        dependencies = cls._collect_dependencies(nodes, {})

        for (filename, nodes) in input_files.iteritems():
            if (filename in output_files):
                producer = output_files[filename][0]
                for consumer in nodes:
                    if producer not in dependencies[consumer]:
                        yield "Node depends on dynamically created file, but not on the node creating it:" + \
                            "\n\tDependent node: %s\n\tFilename: %s\n\tCreated by: %s" \
                            % (consumer, filename, producer)
            elif not os.path.exists(filename):
                yield "Required file does not exist, and is not created by a node:" + \
                            "\n\tFilename: %s\n\tDependent node(s): %s" \
                            % (filename,    "\n\t                   ".join(map(str, nodes)))


    @classmethod
    def _collect_dependencies(cls, nodes, dependencies):
        for node in nodes:
            if node not in dependencies:
                subnodes = node.subnodes | node.dependencies
                if not subnodes:
                    dependencies[node] = frozenset()
                    continue

                cls._collect_dependencies(subnodes, dependencies)
                
                collected = set(subnodes)
                for subnode in subnodes:
                    collected.update(dependencies[subnode])
                dependencies[node] = frozenset(collected)

        return dependencies


    @classmethod
    def _collect_reverse_dependencies(cls, lst, rev_dependencies):
        for node in lst:
            rev_dependencies[node] 
            for dependency in (node.dependencies | node.subnodes):
                rev_dependencies[dependency].add(node)
            cls._collect_reverse_dependencies(node.dependencies, rev_dependencies)
            cls._collect_reverse_dependencies(node.subnodes, rev_dependencies)
        
        

