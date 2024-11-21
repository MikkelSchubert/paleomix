#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import collections
import errno
import fnmatch
import logging
import os
from enum import Enum
from itertools import islice
from shlex import quote
from typing import TYPE_CHECKING, Iterable, TypeVar

from paleomix.common.fileutils import missing_executables
from paleomix.common.versions import Requirement, RequirementError

if TYPE_CHECKING:
    from paleomix.node import Node as Task  # FIXME


class NodeGraphError(RuntimeError):
    pass


class FileStatusCache:
    """Cache used to avoid repeatedly checking the state (existence / mtime) of
    files required / generated by nodes. A new cache is generated for every
    operation (e.g. refreshing all states / manually setting the state of a
    node) to avoid relying on the filesystem staying consistent for long
    periods of time.
    """

    def __init__(self) -> None:
        self._abs_cache: dict[str, str] = {}
        self._mtime_ns_cache: dict[str, int | None] = {}

    def abspath(self, fpath: str) -> str:
        abspath = self._abs_cache.get(fpath)
        if abspath is None:
            abspath = self._abs_cache[fpath] = os.path.abspath(fpath)

        return abspath

    def any_missing_files(self, filepaths: Iterable[str]) -> bool:
        return any(self._mtime_ns(fpath) is None for fpath in filepaths)

    def missing_files(self, filepaths: Iterable[str]) -> Iterable[str]:
        for fpath in filepaths:
            if self._mtime_ns(fpath) is None:
                yield fpath

    def newest_mtime_ns(self, filepaths: Iterable[str]) -> int | None:
        newest = None
        for fpath in filepaths:
            mtime = self._mtime_ns(fpath)
            if mtime is not None and (newest is None or newest < mtime):
                newest = mtime

        return newest

    def oldest_mtime_ns(self, filepaths: Iterable[str]) -> int | None:
        oldest = None
        for fpath in filepaths:
            mtime = self._mtime_ns(fpath)
            if mtime is not None and (oldest is None or oldest > mtime):
                oldest = mtime

        return oldest

    def _mtime_ns(self, fpath: str) -> int | None:
        """Returns the mtime of a path, or None if the path does not exist."""
        if fpath not in self._mtime_ns_cache:
            try:
                mtime = os.stat(fpath).st_mtime_ns
            except OSError as error:
                if error.errno != errno.ENOENT:
                    raise
                mtime = None
            self._mtime_ns_cache[fpath] = mtime
        return self._mtime_ns_cache[fpath]


class StatusEnum(Enum):
    DONE = "done"
    RUNNING = "running"
    RUNNABLE = "runnable"
    QUEUED = "queued"
    ERROR = "error"


class CleanupStrategy(Enum):
    DELETE = "delete"
    KEEP = "keep"
    REQUIRE = "require"

    def __str__(self) -> str:
        # Use `X` rather than `CleanupStrategy.X` when used with argparse, etc.
        return self.value


class TaskStatus:
    obj: Task
    status: StatusEnum
    outdated: bool
    dependencies: set[TaskStatus]
    rev_dependencies: set[TaskStatus]

    # Normalized input/output file paths
    input_files: frozenset[str]
    output_files: frozenset[str]
    auxiliary_files: frozenset[str]
    intermediate_output_files: set[str]

    def __init__(self, task: Task, fscache: FileStatusCache) -> None:
        self.obj = task

        self.status = StatusEnum.DONE
        self.dependencies = set()
        self.rev_dependencies = set()

        self.input_files = frozenset(map(fscache.abspath, task.input_files))
        self.output_files = frozenset(map(fscache.abspath, task.output_files))
        self.auxiliary_files = frozenset(map(fscache.abspath, task.auxiliary_files))
        self.intermediate_output_files = {
            fscache.abspath(filepath) for filepath in task.intermediate_output_files
        }

    def is_queued_and_runnable(self) -> bool:
        if self.status != StatusEnum.QUEUED:
            return False

        return all(task.is_done() for task in self.dependencies)

    def is_needed_by_dependencies(self) -> bool:
        if self.status != StatusEnum.DONE:
            return True

        return any(task.status != StatusEnum.DONE for task in self.rev_dependencies)

    def is_done(self) -> bool:
        return self.status == StatusEnum.DONE

    def __str__(self) -> str:
        return str(self.obj)

    def __repr__(self) -> str:
        return f"TaskStatus({self.obj})"


class NodeGraph:
    tasks: frozenset[Task]
    requirements: tuple[Requirement, ...]

    DONE = StatusEnum.DONE
    RUNNING = StatusEnum.RUNNING
    RUNNABLE = StatusEnum.RUNNABLE
    QUEUED = StatusEnum.QUEUED
    ERROR = StatusEnum.ERROR

    _log: logging.Logger
    _status: dict[Task, TaskStatus]
    _status_counts: dict[StatusEnum, int]
    _intermediate_files: list[str]
    _intermediate_files_strategy: CleanupStrategy

    def __init__(
        self,
        tasks: Iterable[Task],
        fscache: FileStatusCache | None = None,
        intermediate_files: CleanupStrategy = CleanupStrategy.DELETE,
        required_files: Iterable[str] = (),
    ) -> None:
        if not fscache:
            fscache = FileStatusCache()

        self._log = logging.getLogger(__name__)
        self._intermediate_files = []
        self._intermediate_files_strategy = intermediate_files

        self._log.debug("Flattening task graph")
        self._status = self._create_task_status_map(tasks, fscache)
        self.tasks = frozenset(self._status)

        status = self._status.values()
        self._log.debug("Creating producer map linking output files to tasks")
        producers = self._collect_file_producers(status)
        self._log.debug("Linking (reverse) dependent tasks by input/output files")
        self._update_dependencies(status, producers)

        self._log.debug("Determining tasks that must be (re-)run due to outdated files")
        self._flag_outdated_tasks(status, fscache)
        self._log.debug("Mark user-requested intermediate files")
        self._require_intermediate_files(producers, required_files)
        self._log.debug("Iteratively resolving status of tasks")
        self._resolve_task_status(status, producers, fscache)
        self._log.debug("Re-add explicit (validation) dependencies")
        self._link_explicit_dependencies(self._status)
        self._log.debug("Check for intermediate files without dependencies")
        self._check_intermediate_files(self._status)
        self._log.debug("Locating first batch of runnable tasks")
        self._mark_runnable_tasks(self._status)
        self._log.debug("Gathering task state statistics")
        self._status_counts = dict.fromkeys(StatusEnum, 0)
        for it in self._status.values():
            self._status_counts[it.status] += 1

        self._log.debug("Collecting software requirements")
        self.requirements = self._collect_requirements(self.tasks)

    def get_node_states(self) -> dict[Task, TaskStatus]:
        return dict(self._status.items())

    def get_node_state(self, task: Task) -> StatusEnum:
        return self._status[task].status

    def set_node_status(
        self,
        task: Task,
        status: StatusEnum,
    ) -> Iterable[tuple[Task, StatusEnum, StatusEnum]]:
        it = self._status[task]
        old_status = it.status
        it.status = status

        self._status_counts[old_status] -= 1
        self._status_counts[status] += 1

        if status == StatusEnum.RUNNING:
            if old_status != StatusEnum.RUNNABLE:
                raise ValueError(status)  # FIXME
        elif status == StatusEnum.DONE:
            for rdep in it.rev_dependencies:
                if rdep.is_queued_and_runnable():
                    self._log.debug("[%s] %s", StatusEnum.RUNNABLE, rdep)
                    yield rdep.obj, StatusEnum.QUEUED, StatusEnum.RUNNABLE

                    rdep.status = StatusEnum.RUNNABLE
                    self._status_counts[StatusEnum.QUEUED] -= 1
                    self._status_counts[StatusEnum.RUNNABLE] += 1

            # Locate no longer needed output files
            for dependency in it.dependencies:
                self._locate_unneeded_files(dependency)

            self._locate_unneeded_files(it)
        elif status == StatusEnum.ERROR:
            for rdep in it.rev_dependencies:
                yield from self.set_node_status(rdep.obj, status)
        elif status == StatusEnum.RUNNABLE:
            if old_status != StatusEnum.RUNNING:
                raise ValueError(status)  # FIXME
        else:
            raise ValueError(status)

        yield task, old_status, status

    def get_state_counts(self) -> dict[StatusEnum, int]:
        return dict(self._status_counts)

    def get_and_reset_intermediate_files(self) -> list[str]:
        self._intermediate_files, intermediate_files = [], self._intermediate_files

        return intermediate_files

    @classmethod
    def check_version_requirements(
        cls,
        requirements: Iterable[Requirement],
        *,
        force: bool = False,
    ) -> bool:
        executables = {requirement.executable for requirement in requirements}
        missing_execs = missing_executables(executables - {None})

        any_errors = False
        log = logging.getLogger(__name__)
        for requirement in requirements:
            name = requirement.name
            if requirement.executable in missing_execs:
                if requirement.name != requirement.executable:
                    name = f"Executable for {name} ({requirement.executable!r})"

                log.error(" [☓] %s not found", name)
                any_errors = True
                continue

            try:
                version = requirement.version(force=force)
                if version:
                    name = f"{name} v{version}"

                if requirement.check():
                    log.info("  [✓] %s ", name)
                else:
                    any_errors = True
                    log.error(
                        " [☓] %s found but %s is required", name, requirement.specifiers
                    )
            except OSError as error:
                any_errors = True
                log.error(" [☓] %s: %s", name, error)
            except RequirementError as error:
                any_errors = True
                log.error(" [☓] %s: %s", name, "\n     ".join(str(error).split("\n")))

        if any_errors:
            log.error(
                "Please refer to the PALEOMIX installation instructions at\n"
                "  https://paleomix.readthedocs.io/en/stable/"
            )

        return not any_errors

    def check_file_dependencies(self, fscache: FileStatusCache) -> bool:
        in_files: dict[str, list[TaskStatus]] = collections.defaultdict(list)
        aux_files: dict[str, list[TaskStatus]] = collections.defaultdict(list)
        out_files: set[str] = set()

        for task in self._status.values():
            out_files.update(task.output_files)
            for fpath in task.input_files:
                in_files[fpath].append(task)
            for fpath in task.auxiliary_files:
                aux_files[fpath].append(task)

        any_missing_files = False
        for label, tasks, required_files in [
            ("input", in_files, in_files.keys() - out_files),
            ("auxiliary", aux_files, aux_files),
        ]:
            missing_files = sorted(fscache.missing_files(required_files))
            if missing_files:
                any_missing_files = True
                self._log.error("Required %s files not found:", label)
                for fpath in missing_files:
                    self._log.error("  - %s", quote(fpath))

                    for task in _summarize(tasks[fpath]):
                        self._log.debug("    required when %s", task)
                        # Failed tasks are marked, to be able to report initial states
                        for _ in self.set_node_status(task.obj, StatusEnum.ERROR):
                            pass

        return not any_missing_files

    def _create_task_status_map(
        self,
        graph: Iterable[Task],
        fscache: FileStatusCache,
    ) -> dict[Task, TaskStatus]:
        queue = set(graph)
        tasks: dict[Task, TaskStatus] = {}

        while queue:
            task = queue.pop()
            tasks[task] = status = TaskStatus(task=task, fscache=fscache)
            for dependency in task.dependencies:
                if dependency not in tasks:
                    queue.add(dependency)

            if self._intermediate_files_strategy == CleanupStrategy.REQUIRE:
                status.intermediate_output_files.clear()

        return tasks

    def _collect_file_producers(
        self, tasks: Iterable[TaskStatus]
    ) -> dict[str, TaskStatus]:
        """Collects actual dependencies, as determined input files used by each task"""
        producers: dict[str, TaskStatus] = {}
        clobbers: dict[str, set[TaskStatus]] = collections.defaultdict(set)
        for task in tasks:
            for output_file in task.output_files:
                existing_producer = producers.get(output_file)
                if existing_producer is not None:
                    clobbers[output_file].add(existing_producer)
                    clobbers[output_file].add(task)

                producers[output_file] = task

        if clobbers:
            self._on_error_producers_clobber(clobbers)

        return producers

    def _on_error_producers_clobber(
        self,
        clobbers: dict[str, set[TaskStatus]],
        max_items: int = 4,
    ) -> None:
        self._log.error("Multiple tasks write to the same output files:")
        for fpath in _summarize(clobbers, max_items=max_items):
            self._log.error("  File %s", quote(fpath))

            tasks = clobbers[fpath]
            for task in _summarize(tasks, max_items=max_items):
                self._log.error("    Created by %s", task)

            if len(tasks) > max_items:
                self._log.error("      and %s more", len(tasks) - max_items)

        if len(clobbers) > max_items:
            self._log.error("  and %s more", len(clobbers) - max_items)

        raise NodeGraphError("Aborted due to bug!")

    @staticmethod
    def _update_dependencies(
        tasks: Iterable[TaskStatus],
        producers: dict[str, TaskStatus],
    ) -> None:
        """Collects actual dependencies, as determined input files used by each task"""
        for task in tasks:
            for input_file in task.input_files:
                dependency = producers.get(input_file)
                if dependency is not None:
                    task.dependencies.add(dependency)
                    dependency.rev_dependencies.add(task)

    def _flag_outdated_tasks(
        self,
        tasks: Iterable[TaskStatus],
        fscache: FileStatusCache,
    ) -> None:
        """Determine tasks that must be re-run based on upstream files being newer;
        this does not take missing files into consideration, since these are handled
        iteratively elsewhere.
        """
        cache: dict[TaskStatus, int | None] = {}

        def _calculate_implied_output_age(task: TaskStatus) -> int | None:
            if task not in cache:
                timestamps: list[int] = []
                for dependency in task.dependencies:
                    timestamp = _calculate_implied_output_age(dependency)
                    if timestamp is not None:
                        timestamps.append(timestamp)

                for io_files in (task.input_files, task.output_files):
                    timestamp = fscache.newest_mtime_ns(io_files)
                    if timestamp is not None:
                        timestamps.append(timestamp)

                cache[task] = max(timestamps, default=None)

            return cache[task]

        for task in tasks:
            timestamp = fscache.oldest_mtime_ns(task.output_files)
            if timestamp is not None:
                for dependency in task.dependencies:
                    in_timestamp = _calculate_implied_output_age(dependency)
                    if in_timestamp is not None and in_timestamp > timestamp:
                        self._log.debug(" - [%s] %s", StatusEnum.QUEUED, task)
                        task.status = StatusEnum.QUEUED
                        break

    def _require_intermediate_files(
        self,
        producers: dict[str, TaskStatus],
        globs: Iterable[str],
    ) -> None:
        any_errors = False
        for glob in globs:
            nmarked = 0
            tstatus = None
            for filename in fnmatch.filter(producers, glob):
                tstatus = producers[filename]
                if filename in tstatus.intermediate_output_files:
                    tstatus.intermediate_output_files.remove(filename)
                    nmarked += 1

            self._log.debug("Required %i intermediate files matching %r", nmarked, glob)
            if tstatus is None:
                self._log.error("No output files found for --require-files %r", glob)
                any_errors = True

        if any_errors:
            raise NodeGraphError("Files requested via --require-files do not exist")

    def _resolve_task_status(
        self,
        tasks: Iterable[TaskStatus],
        producers: dict[str, TaskStatus],
        fscache: FileStatusCache,
    ) -> None:
        queue = set(tasks)
        while queue:
            task = queue.pop()

            if task.status != StatusEnum.QUEUED:
                required_files = task.output_files - task.intermediate_output_files
                if not fscache.any_missing_files(required_files):
                    continue

                self._log.debug(" [%s] %s", StatusEnum.QUEUED, task)
                task.status = StatusEnum.QUEUED

            # (Re)generate input files as needed
            for missing_file in fscache.missing_files(task.input_files):
                producer = producers.get(missing_file)
                if producer is not None and producer.status != StatusEnum.QUEUED:
                    self._log.debug(" [%s] %s", StatusEnum.QUEUED, producer)
                    producer.status = StatusEnum.QUEUED
                    queue.add(producer)

            # Update downstream tasks not already queued to be re-run
            for dependency in task.rev_dependencies:
                if dependency.status != StatusEnum.QUEUED:
                    self._log.debug(" [%s] %s", StatusEnum.QUEUED, dependency)
                    dependency.status = StatusEnum.QUEUED
                    queue.add(dependency)

    def _link_explicit_dependencies(self, tasks: dict[Task, TaskStatus]) -> None:
        # Explicit dependencies are those not corresponding to use of output files and
        # typically equate to validation steps, and must therefore be added separately
        # from dependencies determined by I/O.
        def _add_dependencies(dependency: TaskStatus, queue: list[TaskStatus]) -> None:
            while queue:
                task: TaskStatus = queue.pop()
                if dependency not in task.dependencies:
                    if task.is_done():
                        # Dependency must be placed lower in the graph to have an effect
                        queue.extend(task.rev_dependencies)
                        continue

                    self._log.debug(" [extra dependency] %r -> %r", dependency, task)
                    task.dependencies.add(dependency)
                    dependency.rev_dependencies.add(task)

        for task, tstatus in tasks.items():
            for dependency in task.dependencies:
                dependency = tasks[dependency]

                if dependency not in tstatus.dependencies:
                    # The original dependee should always depend on the task
                    tstatus.dependencies.add(dependency)
                    dependency.rev_dependencies.add(tstatus)

                    # If a validation step is to be run, then downstream tasks must be
                    # blocked even if the task depending on the validation step is done.
                    if tstatus.is_done() and not dependency.is_done():
                        _add_dependencies(dependency, list(tstatus.rev_dependencies))

    def _check_intermediate_files(self, tasks: dict[Task, TaskStatus]) -> None:
        any_errors = False
        for task, tstatus in tasks.items():
            if (
                tstatus.intermediate_output_files.issuperset(tstatus.output_files)
                and not tstatus.rev_dependencies
            ):
                self._log.error("Intermediate task with no dependencies: %s", task)
                any_errors = True

        if any_errors:
            raise NodeGraphError("Aborted due to bug!")

    def _mark_runnable_tasks(self, tasks: dict[Task, TaskStatus]) -> None:
        for task in tasks.values():
            if task.status == StatusEnum.QUEUED:
                if all(dep.status == StatusEnum.DONE for dep in task.dependencies):
                    self._log.debug(" [%s] %s", StatusEnum.RUNNABLE, task)
                    task.status = StatusEnum.RUNNABLE
            elif task.status == StatusEnum.DONE:
                self._locate_unneeded_files(task)

    @staticmethod
    def _collect_requirements(tasks: Iterable[Task]) -> tuple[Requirement, ...]:
        executables: set[str] = set()
        requirements: set[Requirement] = set()
        for task in tasks:
            executables.update(task.executables)
            requirements.update(task.requirements)

        # Executables used by requirement checks
        requirement_execs = {requirement.executable for requirement in requirements}

        # Create dummy Requirements object for any executables without Requirements
        for executable in executables - requirement_execs:
            requirement = Requirement(executable)
            # Handle the presence of "%(PYTHON)s"
            if requirement.executable not in requirement_execs:
                requirements.add(requirement)

        return tuple(sorted(requirements, key=lambda req: req.name))

    def _locate_unneeded_files(self, task: TaskStatus) -> None:
        if not task.is_done():
            raise AssertionError("called _locate_unneeded_files on unfinished task")

        if (
            self._intermediate_files_strategy == CleanupStrategy.DELETE
            and task.intermediate_output_files
            and not task.is_needed_by_dependencies()
        ):
            cache = FileStatusCache()
            missing_files = cache.missing_files(task.intermediate_output_files)
            removable_files = task.intermediate_output_files.difference(missing_files)
            self._intermediate_files.extend(removable_files)


T = TypeVar("T")


def _summarize(items: Iterable[T], max_items: int = 4) -> Iterable[T]:
    unique_items = {str(it): it for it in items}
    sorted_items = sorted(unique_items.items())

    return (value for _, value in islice(sorted_items, max_items))
