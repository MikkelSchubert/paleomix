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
import logging
import collections

import pypeline.common.versions as versions

from pypeline.node import MetaNode
from pypeline.common.fileutils import \
     reroot_path, \
     missing_executables
from pypeline.common.utilities import \
     safe_coerce_to_frozenset


# Max number of error messages of each type
_MAX_ERROR_MESSAGES = 10


class NodeGraphError(RuntimeError):
    pass


class NodeGraph:
    NUMBER_OF_STATES = 6
    DONE, RUNNING, RUNABLE, QUEUED, OUTDATED, ERROR \
      = range(NUMBER_OF_STATES)

    def __init__(self, nodes):
        self._state_observers = []
        self._states = {}

        nodes = safe_coerce_to_frozenset(nodes)

        self._logger = logging.getLogger(__name__)
        self._reverse_dependencies = collections.defaultdict(set)
        self._collect_reverse_dependencies(nodes, self._reverse_dependencies, set())
        self._intersections = self._calculate_intersections()
        self._top_nodes = [node for (node, rev_deps) in self._reverse_dependencies.iteritems() if not rev_deps]

        self._logger.info("  - Checking file dependencies ...")
        self._check_file_dependencies(self._reverse_dependencies)
        self._logger.info("  - Checking for required executables ...")
        self._check_required_executables(self._reverse_dependencies)
        self._logger.info("  - Checking version requirements ...")
        self._check_version_requirements(self._reverse_dependencies)
        self._logger.info("  - Determining states ...")
        self.refresh_states()
        self._logger.info("  - Ready ...\n")


    def get_node_state(self, node):
        return self._states[node]


    def set_node_state(self, node, state):
        if state not in (NodeGraph.RUNNING, NodeGraph.ERROR, NodeGraph.DONE):
            raise ValueError("Cannot set states other than RUNNING and ERROR, or DONE.")
        old_state = self._states[node]
        if state == old_state:
            return

        self._states[node] = state
        self._notify_state_observers(node, old_state, state, True)

        intersections = dict(self._intersections[node])

        # Not all nodes may need to be updated, but we still need to
        # traverse the "graph" (using the intersection counts) in order
        # to ensure that all nodes that need to be updated are updated.
        requires_update = dict.fromkeys(intersections, False)
        for dependency in self._reverse_dependencies[node]:
            requires_update[dependency] = True

        while any(requires_update.itervalues()):
            for (node, count) in intersections.items():
                if not count:
                    has_changed = False
                    if requires_update[node]:
                        old_state = self._states.pop(node)
                        new_state = self._update_node_state(node)
                        if (new_state != old_state):
                            self._notify_state_observers(node, old_state, new_state, False)
                            has_changed = True

                    for dependency in self._reverse_dependencies[node]:
                        intersections[dependency] -= 1
                        requires_update[dependency] |= has_changed

                    intersections.pop(node)
                    requires_update.pop(node)


    def __iter__(self):
        """Returns a graph of nodes."""
        return iter(self._top_nodes)


    def iterflat(self):
        return iter(self._reverse_dependencies)


    def refresh_states(self):
        states = {}
        for (node, state) in self._states.iteritems():
            if state in (self.ERROR, self.RUNNING):
                states[node] = state
        self._states = states
        for node in self._reverse_dependencies:
            self._update_node_state(node)
        self._refresh_state_observers()


    def add_state_observer(self, observer):
        """Add an observer of changes to the node-graph. The observer
        is expected to have the following functions:

        refresh(nodegraph):
          Called when an observer has been added, or when 'refresh_states'
          has been called on the nodegraph. The observer should rebuild any
          internal state at this point.

        state_changed(node, old_state, new_state, is_primary):
          Called when the state of an node has changed. 'is_primary' is
          True only if the node for which 'set_node_state' was called,
          and false for nodes the state of which changed as a consequence
          of the change to the node marked 'is_primary'. This includes
          ERROR propegating, MetaNodes being DONE when all their subnodes
          are DONE and more."""
        self._state_observers.append(observer)
        observer.refresh(self)


    def _notify_state_observers(self, node, old_state, new_state, is_primary):
        for observer in self._state_observers:
            observer.state_changed(node, old_state, new_state, is_primary)


    def _refresh_state_observers(self):
        for observer in self._state_observers:
            observer.refresh(self)


    def _calculate_intersections(self):
        def count_nodes(node, counts):
            for node in self._reverse_dependencies[node]:
                if node in counts:
                    counts[node] += 1
                else:
                    counts[node] = 1
                    count_nodes(node, counts)
            return counts

        intersections = {}
        for node in self._reverse_dependencies:
            counts = count_nodes(node, {})
            for dependency in self._reverse_dependencies[node]:
                counts[dependency] -= 1
            intersections[node] = counts

        return intersections


    def _update_node_state(self, node):
        if node in self._states:
            return self._states[node]

        # Update sub-nodes before checking for fixed states
        dependency_states = set((NodeGraph.DONE,))
        for dependency in node.dependencies:
            dependency_states.add(self._update_node_state(dependency))

        subnode_states = set()
        for subnode in node.subnodes:
            subnode_states.add(self._update_node_state(subnode))

        try:
            state = max(subnode_states | dependency_states)
            if isinstance(node, MetaNode):
                if (state < NodeGraph.ERROR):
                    # Mark a meta-node as running if a subnode is running
                    running_states = (NodeGraph.RUNNING in subnode_states)
                    if (state == NodeGraph.RUNABLE):
                        # MetaNodes cannot actually be run!
                        state = NodeGraph.QUEUED
                    if (state == NodeGraph.RUNNING) and not running_states:
                        # A dependency is running, don't inherit state
                        state = NodeGraph.QUEUED
                    elif running_states:
                        # Prevent RUNABLE/QUEUED/OUTDATED from taking precedence
                        state = NodeGraph.RUNNING
            elif state == NodeGraph.DONE:
                if not node.is_done or node.is_outdated:
                    state = NodeGraph.RUNABLE
            elif state in (NodeGraph.RUNNING, NodeGraph.RUNABLE, NodeGraph.QUEUED):
                if node.is_done:
                    state = NodeGraph.OUTDATED
                else:
                    state = NodeGraph.QUEUED
        except OSError, error:
            # Typically hapens if base input files are removed, causing a node that
            # 'is_done' to call modified_after on missing files in 'is_outdated'
            self._logger.error("OSError checking state of Node: %s", error)
            state = NodeGraph.ERROR
        self._states[node] = state

        return state


    @classmethod
    def _check_required_executables(cls, nodes):
        exec_filenames = set()
        for node in nodes:
            exec_filenames.update(node.executables)

        missing_exec   = missing_executables(exec_filenames)
        if missing_exec:
            raise NodeGraphError("Required executables are missing:\n\t%s" \
                                % ("\n\t".join(sorted(missing_exec))))


    @classmethod
    def _check_version_requirements(cls, nodes):
        exec_requirements = set()
        for node in nodes:
            exec_requirements.update(node.requirements)

        try:
            for requirement in exec_requirements:
                requirement()
        except versions.VersionRequirementError, error:
            raise NodeGraphError(error)
        except OSError, error:
            raise NodeGraphError("Could not check version requirements for %s:\n\t%s" \
                                 % (requirement.name, error))

    @classmethod
    def _check_file_dependencies(cls, nodes):
        files = ("input_files", "output_files")
        files = dict((key, collections.defaultdict(set)) for key in files)
        # Auxiliary files are treated as input files
        files["auxiliary_files"] = files["input_files"]

        for node in nodes:
            for (attr, nodes_by_file) in files.iteritems():
                for filename in getattr(node, attr):
                    nodes_by_file[filename].add(node)

        max_messages = range(_MAX_ERROR_MESSAGES)
        error_messages = []
        error_messages.extend(zip(max_messages, cls._check_output_files(files["output_files"])))
        error_messages.extend(zip(max_messages, cls._check_input_dependencies(files["input_files"],
                                                                              files["output_files"], nodes)))

        if error_messages:
            messages = []
            for (_, error) in error_messages:
                for line in error.split("\n"):
                    messages.append("\t" + line)

            raise NodeGraphError("Errors detected during graph construction (max %i shown):\n%s" \
                                % (_MAX_ERROR_MESSAGES * 2, "\n".join(messages)),)


    @classmethod
    def _check_output_files(cls, output_files):
        """Checks dict of output files to nodes for cases where
        multiple nodes create the same output file.

        The directory component of paths are realized in order to
        detect cases where nodes create the same file, but via
        different paths (e.g. due to relative/absolute paths, or
        due to use of symbolic links). Since output files are
        replaced, not modified in place, it is not nessesary to
        compare files themselves."""
        dirpath_cache, real_output_files = {}, {}
        for (filename, nodes) in output_files.iteritems():
            dirpath = os.path.dirname(filename)
            if dirpath not in dirpath_cache:
                dirpath_cache[dirpath] = os.path.realpath(dirpath)

            real_output_file = reroot_path(dirpath_cache[dirpath], filename)
            real_output_files.setdefault(real_output_file, []).extend(nodes)

        for (filename, nodes) in real_output_files.iteritems():
            if (len(nodes) > 1):
                nodes = _summarize_nodes(nodes)
                yield "Multiple nodes create the same (clobber) output-file:" \
                      "\n\tFilename: %s\n\tNodes: %s" \
                      % (filename, "\n\t       ".join(nodes))


    @classmethod
    def _check_input_dependencies(cls, input_files, output_files, nodes):
        dependencies = cls._collect_dependencies(nodes, {})

        for (filename, nodes) in sorted(input_files.items(), key = lambda v: v[0]):
            if (filename in output_files):
                producers = output_files[filename]
                bad_nodes = set()
                for consumer in nodes:
                    if not (producers & dependencies[consumer]):
                        bad_nodes.add(consumer)

                if bad_nodes:
                    producer  = iter(producers).next()
                    bad_nodes = _summarize_nodes(bad_nodes)
                    yield "Node depends on dynamically created file, but not on the node creating it:" + \
                                "\n\tFilename: %s\n\tCreated by: %s\n\tDependent node(s): %s" \
                                % (filename, producer, "\n\t                   ".join(bad_nodes))
            elif not os.path.exists(filename):
                nodes = _summarize_nodes(nodes)
                yield "Required file does not exist, and is not created by a node:" + \
                            "\n\tFilename: %s\n\tDependent node(s): %s" \
                            % (filename,    "\n\t                   ".join(nodes))


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
    def _collect_reverse_dependencies(cls, lst, rev_dependencies, processed):
        for node in lst:
            if node not in processed:
                processed.add(node)

                # Initialize default-dict
                rev_dependencies[node] # pylint: disable=W0104

                subnodes = node.dependencies | node.subnodes
                for dependency in subnodes:
                    rev_dependencies[dependency].add(node)
                cls._collect_reverse_dependencies(subnodes, rev_dependencies, processed)



def _summarize_nodes(nodes):
    nodes = list(sorted(set(map(str, nodes))))
    if len(nodes) > 4:
        nodes = nodes[:5] + ["and %i more nodes ..." % len(nodes)]
    return nodes
