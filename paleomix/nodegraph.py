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
import os

import paleomix.common.versions as versions

from paleomix.common.fileutils import missing_executables, missing_files
from paleomix.common.utilities import safe_coerce_to_frozenset


class FileStatusCache:
    """Cache used to avoid repeatedly checking the state (existance / mtime) of
    files required / generated by nodes. A new cache is generated for every
    operation (e.g. refreshing all states / manually setting the state of a
    node) to avoid relying on the filesystem staying consistant for long
    periods of time.
    """

    def __init__(self):
        self._stat_cache = {}

    def files_exist(self, fpaths):
        """Returns true if all paths listed in fpaths exist."""
        return all((self._get_state(fpath) is not None) for fpath in fpaths)

    def missing_files(self, fpaths):
        """Returns a list of paths in fpaths that do not exist."""
        return [fpath for fpath in fpaths if (self._get_state(fpath) is None)]

    def are_files_outdated(self, input_files, output_files):
        """Returns true if any 'input' files have a time-stamp that post-date
        any time-stamp for the 'output' files, indicating that one or more of
        the 'input' files have changed since the creation of the 'output'.

        The function also returns true if any files are missing, as this would
        indicate either the 'output' files or both the 'input' and 'output'
        files would need to be rebuilt.
        """
        input_timestamps = []
        if not self._get_states(input_files, input_timestamps):
            return True

        output_timestamps = []
        if not self._get_states(output_files, output_timestamps):
            return True

        return max(input_timestamps) > min(output_timestamps)

    def _get_states(self, filenames, dst):
        """Collects the mtimes for a set of filenames, returning true if all
        could be collected, and aborting early and returning false otherwise.
        """
        for filename in filenames:
            timestamp = self._get_state(filename)
            if timestamp is None:
                return False

            dst.append(timestamp)

        return True

    def _get_state(self, fpath):
        """Returns the mtime of a path, or None if the path does not exist."""
        if fpath not in self._stat_cache:
            try:
                mtime = os.path.getmtime(fpath)
            except OSError as error:
                if error.errno != errno.ENOENT:
                    raise
                mtime = None
            self._stat_cache[fpath] = mtime
        return self._stat_cache[fpath]


class NodeGraphError(RuntimeError):
    pass


class NodeGraph:
    NUMBER_OF_STATES = 6
    DONE, RUNNING, RUNABLE, QUEUED, OUTDATED, ERROR = range(NUMBER_OF_STATES)

    def __init__(
        self,
        nodes,
        allow_missing_files=False,
        cache_factory=FileStatusCache,
    ):
        self._cache_factory = cache_factory
        self._states = {}
        self._state_counts = [0] * self.NUMBER_OF_STATES

        nodes = safe_coerce_to_frozenset(nodes)

        self._logger = logging.getLogger(__name__)
        self._reverse_dependencies = collections.defaultdict(set)
        self._collect_reverse_dependencies(nodes, self._reverse_dependencies, set())
        self._intersections = {}
        self._top_nodes = [
            node
            for (node, rev_deps) in self._reverse_dependencies.items()
            if not rev_deps
        ]

        self._logger.info("Checking file dependencies")
        if not self._check_file_dependencies(self._reverse_dependencies):
            raise NodeGraphError("Aborting due to input/output file error")

        self._logger.info("Checking for auxiliary files")
        if not self._check_auxiliary_files(self._reverse_dependencies):
            if not allow_missing_files:
                raise NodeGraphError(
                    "Please refer to the PALEOMIX installation instructions at "
                    "https://paleomix.readthedocs.io/en/stable/"
                )

        self._logger.info("Checking required software")
        if not self._check_version_requirements(self._reverse_dependencies):
            if not allow_missing_files:
                raise NodeGraphError(
                    "Please refer to the PALEOMIX installation instructions at "
                    "https://paleomix.readthedocs.io/en/stable/"
                )

        self._logger.info("Determining states")
        self._refresh_states()
        self._logger.info("Ready")

    def get_node_state(self, node):
        return self._states[node]

    def set_node_state(self, node, state):
        if state not in (NodeGraph.RUNNING, NodeGraph.ERROR, NodeGraph.DONE):
            raise ValueError("Invalid state: %r" % (state,))
        old_state = self._states[node]
        if state == old_state:
            return

        self._states[node] = state
        self._state_counts[old_state] -= 1
        self._state_counts[state] += 1

        self._log_node_changes(node, old_state, state)

        intersections = self._calculate_intersections(node)

        # Not all nodes may need to be updated, but we still need to
        # traverse the "graph" (using the intersection counts) in order
        # to ensure that all nodes that need to be updated are updated.
        requires_update = dict.fromkeys(intersections, False)
        for dependency in self._reverse_dependencies[node]:
            requires_update[dependency] = True

        cache = self._cache_factory()
        while any(requires_update.values()):
            for (node, count) in tuple(intersections.items()):
                if not count:
                    has_changed = False
                    if requires_update[node]:
                        old_state = self._states.pop(node)
                        new_state = self._update_node_state(node, cache)
                        has_changed |= new_state != old_state

                        self._state_counts[old_state] -= 1
                        self._state_counts[new_state] += 1

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

    def get_state_counts(self):
        return list(self._state_counts)

    def _refresh_states(self):
        states = {}
        cache = self._cache_factory()
        for (node, state) in self._states.items():
            if state in (self.ERROR, self.RUNNING):
                states[node] = state
        self._states = states

        for node in self._reverse_dependencies:
            self._update_node_state(node, cache)

        state_counts = [0] * self.NUMBER_OF_STATES
        for state in states.values():
            state_counts[state] += 1
        self._state_counts = state_counts

    def _log_node_changes(self, node, old_state, new_state):
        if new_state in (self.RUNNING, self.DONE):
            running = self._state_counts[self.RUNNING]
            remaining = (
                sum(self._state_counts)
                - self._state_counts[self.DONE]
                - self._state_counts[self.ERROR]
            )

            if new_state == self.RUNNING:
                self._logger.info("[%i/%i] Started %s", running, remaining, node)
            elif new_state == self.DONE:
                self._logger.info("[%i/%i] Finished %s", running, remaining, node)

    def _calculate_intersections(self, for_node):
        def count_nodes(node, counts):
            for node in self._reverse_dependencies[node]:
                if node in counts:
                    counts[node] += 1
                else:
                    counts[node] = 1
                    count_nodes(node, counts)
            return counts

        if for_node not in self._intersections:
            counts = count_nodes(for_node, {})
            for dependency in self._reverse_dependencies[for_node]:
                counts[dependency] -= 1
            self._intersections[for_node] = counts

        return dict(self._intersections[for_node])

    def _update_node_state(self, node, cache):
        if node in self._states:
            return self._states[node]

        # Update sub-nodes before checking for fixed states
        dependency_states = set((NodeGraph.DONE,))
        for dependency in node.dependencies:
            dependency_states.add(self._update_node_state(dependency, cache))

        state = max(dependency_states)
        if state == NodeGraph.DONE:
            if not self.is_done(node, cache):
                state = NodeGraph.RUNABLE
            elif not cache.files_exist(node.input_files):
                # Somehow the input files have gone missing, despite the
                # dependant nodes being done; this implies this node is
                # outdated, since the input-files should be re-generated, but
                # obviously it is not possible to run it at this point.
                missing = cache.missing_files(node.input_files)
                self._logger.error(
                    "Input file(s) missing for node; may have been moved while the "
                    "pipeline was running. Cannot proceed:\n"
                    "    Node = %s\n    Files = %s\n",
                    node,
                    "\n            ".join(missing),
                )
                state = NodeGraph.ERROR
            elif self.is_outdated(node, cache):
                state = NodeGraph.RUNABLE
        elif state in (NodeGraph.RUNNING, NodeGraph.RUNABLE, NodeGraph.QUEUED):
            if self.is_done(node, cache):
                state = NodeGraph.OUTDATED
            else:
                state = NodeGraph.QUEUED
        self._states[node] = state

        return state

    @classmethod
    def is_done(cls, node, cache):
        """Returns true if the node itself is done; this only implies that the
        output files generated by this node exists. The files themselves may
        be outdated.
        """
        return cache.files_exist(node.output_files)

    @classmethod
    def is_outdated(cls, node, cache):
        """Returns true if the not is not done or if one or more of the input
        files appear to have been changed since the creation of the output
        files (based on the timestamps). A node that lacks either input or
        output files is never considered outdated.
        """
        if not (node.input_files and node.output_files):
            return False

        return cache.are_files_outdated(node.input_files, node.output_files)

    def _check_required_executables(self, nodes):
        exec_filenames = set()
        for node in nodes:
            exec_filenames.update(node.executables)

            # Requirements may include executables not invoked directly
            for requirement in node.requirements:
                if isinstance(requirement, versions.RequirementObj):
                    executable = requirement.executable
                    if executable is not None:
                        exec_filenames.add(executable)

        missing_exec = missing_executables(exec_filenames)
        if missing_exec:
            self._logger.error("Required executable(s) not found:")
            for name in sorted(missing_exec):
                self._logger.error(" - %s", name)

        return not missing_exec

    def _check_version_requirements(self, nodes):
        if not self._check_required_executables(nodes):
            return False

        exec_requirements = set()
        for node in nodes:
            exec_requirements.update(node.requirements)

        def _key_func(reqobj):
            # Sort priority in decreasing order, name in increasing order
            return (-reqobj.priority, reqobj.name)

        any_errors = False
        for requirement in sorted(exec_requirements, key=_key_func):
            try:
                name = requirement.name
                version = ".".join(str(value) for value in requirement.version)
                if version:
                    name = "%s v%s" % (name, version)

                self._logger.info(" - Found %s", name)

                requirement()
            except versions.VersionRequirementError as error:
                any_errors = True
                self._logger.error(error)
            except OSError as error:
                any_errors = True
                self._logger.error(
                    "Could not check version for %s:\n\t%s" % (requirement.name, error)
                )

        return not any_errors

    def _check_file_dependencies(self, nodes, max_errors=10):
        input_files = collections.defaultdict(set)
        output_files = collections.defaultdict(set)
        for node in nodes:
            for filename in node.input_files:
                input_files[filename].add(node)

            for filename in node.output_files:
                output_files[filename].add(node)

        any_errors = False
        if not self._check_output_files(output_files, max_errors):
            any_errors = True

        if not self._check_input_files(input_files, output_files, nodes, max_errors):
            any_errors = True

        return not any_errors

    def _check_output_files(self, output_files, max_errors=10):
        """Checks dict of output files to nodes for cases where
        multiple nodes create the same output file.

        The directory component of paths are realized in order to
        detect cases where nodes create the same file, but via
        different paths (e.g. due to relative/absolute paths, or
        due to use of symbolic links). Since output files are
        replaced, not modified in place, it is not nessesary to
        compare files themselves."""
        dirpath_cache = {}
        filepaths = collections.defaultdict(list)
        for filename, nodes in output_files.items():
            dirpath, filename = os.path.split(filename)
            if dirpath not in dirpath_cache:
                dirpath_cache[dirpath] = os.path.realpath(dirpath)

            filepath = os.path.join(dirpath_cache[dirpath], filename)
            filepaths[filepath].extend(nodes)

        clobbered_files = {
            key: values for key, values in filepaths.items() if len(values) > 1
        }

        if clobbered_files:
            for filename, nodes in sorted(clobbered_files.items()):
                self._logger.error("Multiple nodes write to file %r:", filename)
                self._logger.debug("depended on by ")
                for line in _summarize_nodes(nodes):
                    self._logger.debug("  %s", line)

                self._logger.error(
                    "This is probably a bug in PALEOMIX; please report at "
                    "https://github.com/MikkelSchubert/paleomix/issues/new"
                )

        return not clobbered_files

    def _check_input_files(self, input_files, output_files, nodes, max_errors=10):
        dependencies = self._collect_dependencies(nodes, {})
        any_errors = False

        for (filename, nodes) in sorted(input_files.items(), key=lambda v: v[0]):
            if filename in output_files:
                producers = output_files[filename]
                bad_nodes = set()
                for consumer in nodes:
                    if not (producers & dependencies[consumer]):
                        bad_nodes.add(consumer)

                if bad_nodes:
                    any_errors = True
                    producer = next(iter(producers))
                    bad_nodes = _summarize_nodes(bad_nodes)
                    self._logger.error(
                        "Tasks depends on file, but not on the task creating it: %r",
                        filename,
                    )

                    self._logger.debug("  produced by %s", producer)
                    self._logger.debug("  depended on by")
                    for line in _summarize_nodes(bad_nodes):
                        self._logger.debug("    %s", line)

                    self._logger.error(
                        "This is probably a bug in PALEOMIX; please report at "
                        "https://github.com/MikkelSchubert/paleomix/issues/new"
                    )

            elif not os.path.exists(filename):
                any_errors = True
                self._logger.error("Required input file does not exist: %r", filename)
                for line in _summarize_nodes(nodes):
                    self._logger.debug("  required by %s", line)

        return not any_errors

    def _check_auxiliary_files(self, nodes):
        auxiliary_files = set()
        for node in nodes:
            auxiliary_files.update(node.auxiliary_files)

        missing_aux_files = missing_files(auxiliary_files)
        if missing_aux_files:
            self._logger.error("Required files not found:")
            for name in sorted(missing_aux_files):
                self._logger.error(" - %s", name)

        return not missing_aux_files

    @classmethod
    def _collect_dependencies(cls, nodes, dependencies):
        for node in nodes:
            if node not in dependencies:
                subnodes = node.dependencies
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
                rev_dependencies[node]

                subnodes = node.dependencies
                for dependency in subnodes:
                    rev_dependencies[dependency].add(node)
                cls._collect_reverse_dependencies(subnodes, rev_dependencies, processed)


def _summarize_nodes(nodes):
    nodes = list(sorted(set(map(str, nodes))))
    if len(nodes) > 4:
        nodes = nodes[:5] + ["and %i more nodes" % len(nodes)]
    return nodes
