# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import codecs
import contextlib
import json
import logging
import multiprocessing
import os
import os.path
import queue
import signal
import socket
import sys
import time
import traceback
import uuid
from collections import defaultdict
from collections.abc import Collection, Iterable, Iterator, MutableMapping
from multiprocessing import ProcessError
from multiprocessing.connection import Client, Connection, wait
from typing import TYPE_CHECKING, Any, Callable, NoReturn, Optional, Union

import setproctitle
from typing_extensions import Self, TypeAlias

import paleomix
from paleomix.common.procs import (
    RegisteredProcess,
    terminate_all_processes,
    terminate_processes,
)
from paleomix.common.versions import Requirement
from paleomix.core.input import CommandLine, ListTasksEvent, ThreadsEvent
from paleomix.node import Node, NodeError
from paleomix.nodegraph import NodeGraph

EventType: TypeAlias = dict[str, Any]
MessageType: TypeAlias = tuple[int, Optional[BaseException], Optional[list[str]]]
WorkerType: TypeAlias = Union["LocalWorker", "RemoteWorker"]
HandleType: TypeAlias = Union[Connection, socket.socket, int]
QueueType: TypeAlias = "multiprocessing.Queue[MessageType]"

if TYPE_CHECKING:
    AdapterType = logging.LoggerAdapter[logging.Logger]
else:
    AdapterType = logging.LoggerAdapter


# Default location of auto-registration files
AUTO_REGISTER_DIR = os.path.expanduser("~/.paleomix/remote/")


# Protocol overview:
# A.  Connecting to workers:
# A1. Upon connecting, the manager sends a handshake including the working directory,
#     the temporary directory, and any requirements that must be met.
EVT_HANDSHAKE = "HANDSHAKE"
# A2. The worker responds with a handshake. In case of failures, this includes an
#     exception describing that failure, causing the manager to close the connection.
#     The handshake answer is typically followed by a CAPACITY message.
EVT_HANDSHAKE_RESPONSE = "HANDSHAKE_RESPONSE"

# B. Main loop
# B1. Sent to indicate the maximum capacity (currently only threads) has changed.
EVT_CAPACITY = "CAPACITY"
# B2. Sent by the manager to run a task on a worker.
EVT_TASK_START = "TASK_START"
# B3. This message is sent by the worker whenever a task has finished running.
EVT_TASK_DONE = "TASK_DONE"

# C.  Disconnecting:
# C1. Sent by the worker or manager, this message indicates that the worker/manager is
#     shutting down and that any cleanup should be performed before disconnecting.
EVT_SHUTDOWN = "SHUTDOWN"


# Worker (connection) states
_UNINITIALIZED = "uninitialized"
_CONNECTING = "connecting"
_RUNNING = "running"
_TERMINATED = "terminated"


def address_to_name(value: str | tuple[str, int]) -> str:
    if isinstance(value, str):
        return value

    host, port = value

    with contextlib.suppress(OSError):
        host, _ = socket.getnameinfo((host, port), 0)

    return f"{host}:{port}"


class WorkerError(RuntimeError):
    pass


class Manager:
    def __init__(
        self,
        *,
        group: str | None = None,
        threads: int,
        requirements: Collection[Requirement],
        temp_root: str,
        auto_connect: bool = True,
    ) -> None:
        self._local: LocalWorker | None = None
        self._interface = CommandLine()
        self._group = group
        self._workers: dict[str, LocalWorker | RemoteWorker] = {}
        self._json_blacklist: set[str] = set()
        self._worker_blacklist: set[str] = set()
        self._requirements = requirements
        self._threads = threads
        self._temp_root = temp_root
        self._next_auto_connect = 0 if auto_connect else float("inf")
        self._log = logging.getLogger(__name__)

    @property
    def workers(self) -> dict[str, LocalWorker | RemoteWorker]:
        return dict(self._workers)

    @property
    def tasks(self) -> Iterator[Node]:
        for worker in self._workers.values():
            yield from worker.tasks

    def start(self) -> bool:
        if self._local is not None:
            raise WorkerError("Manager already started")

        local_worker = LocalWorker(self, self._threads)
        if not local_worker.connect(self._requirements):
            return False

        self._local = local_worker
        self._workers[self._local.id] = local_worker

        if self._threads <= 0:
            self._log.warning(
                "Local worker process has no threads assigned; either "
                "increase allocation with '+' or start worker processes."
            )

        return self._auto_connect_to_workers()

    def wait_for_workers(self) -> bool:
        self._check_started()
        while not all(worker.status == _RUNNING for worker in self._workers.values()):
            for event in self.poll():
                if (
                    event["event"] == EVT_HANDSHAKE_RESPONSE
                    and event["error"] is not None
                ):
                    return False

        return True

    def shutdown(self) -> None:
        for worker in self._workers.values():
            worker.shutdown()
        self._workers.clear()
        self._local = None

    def poll(self) -> Iterator[EventType]:
        self._check_started()
        self._auto_connect_to_workers()

        for worker, event in self._collect_events():
            if not (isinstance(event, dict) and "event" in event):
                self._log.error("Malformed event from %s: %r", worker.name, event)
                continue

            event["worker"] = worker.id
            event["worker_name"] = worker.name

            if event["event"] == EVT_HANDSHAKE_RESPONSE:
                if event["error"] is not None:
                    self._worker_blacklist.add(worker.id)
                    self._workers.pop(worker.id)
                    self._log.error(
                        "Handshake with worker %r failed:\n  %s",
                        worker.name,
                        event["error"],
                    )

                continue
            elif event["event"] == EVT_SHUTDOWN:
                self._workers.pop(worker.id)

            yield event

        # For now the manager takes care of signalling that there is available capacity,
        # but the plan is for the workers to manage this themselves, so that multiple
        # managers can be served by the same worker
        max_threads = max((it.threads for it in self._workers.values()), default=0)
        for worker in self._workers.values():
            idle_threads = worker.threads - sum(task.threads for task in worker.tasks)
            if idle_threads > 0:
                # Signal that there are available threads
                yield {
                    "event": EVT_CAPACITY,
                    "threads": idle_threads,
                    # The biggest workers are allowed to exceed capacity if needed
                    "overcommit": max_threads and idle_threads == max_threads,
                    "worker": worker.id,
                }

    def start_task(self, worker_id: str, task: Node) -> bool:
        self._check_started()

        worker = self._workers.get(worker_id)
        if worker is None:
            self._log.error("Tried to start task on unknown worker %r", worker_id)
            return False

        return worker.start_task(task, self._temp_root)

    def _collect_events(self) -> list[tuple[WorkerType, EventType]]:
        workers_with_events: dict[WorkerType, list[HandleType]] = defaultdict(list)
        handles: dict[HandleType, WorkerType | None] = {}
        for worker in self._workers.values():
            if worker.has_events():
                self._log.debug("Worker has pending events: %r", worker)
                workers_with_events[worker] = []

            for handle in worker.handles:
                handles[handle] = worker

        # Detect keyboard input from the user (optional)
        for handle in self._interface.handles:
            handles[handle] = None

        # No need to block if we already have events to act on
        timeout = 0 if workers_with_events else 5.0
        for handle in wait(handles, timeout):
            worker = handles[handle]
            self._log.debug("Worker has live handle: %r / %r", worker, handle)
            if worker is None:
                self._handle_keyboard_input()
            else:
                workers_with_events[worker].append(handle)

        events: list[tuple[WorkerType, EventType]] = []
        for worker, live_handles in workers_with_events.items():
            for event in worker.get_events(live_handles):
                self._log.debug("Worker has event: %r / %r", worker, event)
                events.append((worker, event))

        return events

    def _handle_keyboard_input(self) -> None:
        for event in self._interface.process_key_presses():
            if isinstance(event, ThreadsEvent):
                self._handle_threads_event(event.change)
            elif isinstance(event, ListTasksEvent):
                self._handle_list_tasks_event()

    def _handle_threads_event(self, change: int) -> None:
        assert self._local is not None

        threads = min(multiprocessing.cpu_count(), max(0, self._threads + change))
        if threads != self._threads:
            self._log.info("Max threads changed from %i to %i", self._threads, threads)
            self._local.threads = threads
            self._threads = threads

    def _handle_list_tasks_event(self) -> None:
        total_threads = 0
        total_workers = 0
        for worker in self._workers.values():
            tasks = sorted(worker.tasks, key=lambda it: it.id)
            threads = sum(task.threads for task in tasks)

            if tasks:
                self._log.info(
                    "Running %i tasks on %s (using %i/%i threads):",
                    len(tasks),
                    worker.name,
                    threads,
                    worker.threads,
                )

                for idx, task in enumerate(tasks, start=1):
                    self._log.info("  % 2i. %s", idx, task)
            else:
                self._log.info(
                    "No tasks running on %s (using 0/%i threads)",
                    worker.name,
                    worker.threads,
                )

            total_threads += threads
            total_workers += 1

        if total_workers > 1:
            self._log.info(
                "A total of %i threads are used across %i workers",
                total_threads,
                total_workers,
            )

    def _auto_connect_to_workers(self, interval: float = 15.0) -> bool:
        any_errors = False
        if self._next_auto_connect <= time.monotonic():
            for info in self._collect_workers():
                if info["id"] in self._worker_blacklist:
                    continue
                elif info["group"] not in (self._group, None):
                    self._worker_blacklist.add(info["id"])
                    continue

                self._log.info("Connecting to %s:%s", info["host"], info["port"])
                worker = RemoteWorker(**info)
                if worker.id in self._workers:
                    self._log.error("Already connected to worker with id", worker.id)
                    any_errors = True
                elif worker.connect(self._requirements):
                    self._workers[worker.id] = worker
                    # Prevent other pipelines from attempting to connect
                    os.unlink(info["filename"])
                else:
                    self._worker_blacklist.add(info["id"])
                    any_errors = True

            self._next_auto_connect = time.monotonic() + interval

        return not any_errors

    def _collect_workers(self, root: str = AUTO_REGISTER_DIR) -> Iterator[EventType]:
        if os.path.isdir(root):
            for filename in os.listdir(root):
                filename = os.path.join(root, filename)
                if filename in self._json_blacklist or not (
                    filename.endswith(".json") and os.path.isfile(filename)
                ):
                    continue

                try:
                    with open(filename) as handle:
                        data = json.load(handle)

                    data["filename"] = filename
                    data["secret"] = codecs.decode(
                        data["secret"].encode("utf-8"), "base64"
                    )
                except (OSError, KeyError, ValueError) as error:
                    self._log.error("Error reading worker file %r: %s", filename, error)
                    self._json_blacklist.add(filename)
                    continue

                yield data

    def _check_started(self) -> None:
        if self._local is None:
            raise WorkerError("Manager not started")

    def __enter__(self) -> Self:
        self._interface.setup()
        return self

    def __exit__(self, typ: object, exc: object, tb: object) -> None:
        self._interface.teardown()
        self.shutdown()


class RemoteAdapter(AdapterType):
    def __init__(self, name: str, logger: logging.Logger) -> None:
        super().__init__(logger, {"remote": name})

    def process(
        self, msg: logging.LogRecord, kwargs: MutableMapping[str, object]
    ) -> tuple[str, MutableMapping[str, object]]:
        if self.extra is not None:
            return "[{}] {}".format(self.extra["remote"], msg), kwargs

        return str(msg), kwargs


class LocalWorker:
    def __init__(self, manager: Manager, threads: int) -> None:
        self.id = str(uuid.uuid4())
        self._events: list[EventType] = []
        self._manager = manager
        self._queue: QueueType = multiprocessing.Queue()
        self._handles: dict[HandleType, tuple[multiprocessing.Process, Node]] = {}
        self._pending_results: dict[int, Node] = {}
        self._status = _UNINITIALIZED
        self._log = logging.getLogger(__name__)
        self.name = "localhost"
        self.threads = threads

    @property
    def tasks(self) -> Iterator[Node]:
        for _, task in self._handles.values():
            yield task

    @property
    def handles(self) -> Iterator[HandleType]:
        yield from self._handles

    @property
    def status(self) -> str:
        return self._status

    def connect(self, requirements: Collection[Requirement]) -> bool:
        if self._status != _UNINITIALIZED:
            raise WorkerError("Attempted to start already initialized LocalWorker")

        self._log.info("Checking required software on localhost")
        if not NodeGraph.check_version_requirements(requirements):
            return False

        self._status = _RUNNING
        self._events.append({"event": EVT_HANDSHAKE_RESPONSE, "error": None})

        return True

    def start_task(self, task: Node, temp_root: str) -> bool:
        self._check_running()

        self._log.debug("Starting local task %s with id %s", task, task.id)
        proc = RegisteredProcess(
            target=task_wrapper,
            args=(self._queue, task, temp_root),
            daemon=True,
        )

        proc.start()
        self._handles[proc.sentinel] = (proc, task)
        self._pending_results[task.id] = task

        return True

    def has_events(self) -> bool:
        try:
            return bool(self._events) or not self._queue.empty()
        except OSError:
            return False

    def get_events(self, handles: list[HandleType]) -> list[EventType]:
        self._check_running()

        # Messages are given priority, since a process may fail after the message
        events, self._events = self._events, []

        while True:
            try:
                key, error, backtrace = self._queue.get(block=False)
            except queue.Empty:
                break

            # task may have been handled previous, due to the process returning first
            task = self._pending_results.pop(key, None)
            if task is not None:
                events.append(
                    {
                        "event": EVT_TASK_DONE,
                        "task": task,
                        "error": error,
                        "backtrace": backtrace,
                    },
                )

        # clean up wrapper processes
        for handle in handles:
            proc, task = self._handles.pop(handle)

            # `join` should not be required, but is recommended:
            #  https://docs.python.org/3.12/library/multiprocessing.html#all-start-methods
            proc.join()
            if proc.exitcode:
                self._log.debug("Joining local task failed: %s", task)

                # Signal that the task failed, if this has not been done previously
                if self._pending_results.pop(task.id, None):
                    events.append(
                        {
                            "event": EVT_TASK_DONE,
                            "task": task,
                            "error": NodeError(
                                "Worker process for task terminated with exit code "
                                f"{proc.exitcode}: {task}"
                            ),
                            "backtrace": None,
                        },
                    )
            else:
                self._log.debug("Joined local task %s", task)

        return events

    def shutdown(self) -> None:
        if self._status != _TERMINATED:
            self._status = _TERMINATED
            self._log.debug("Shutting down local worker")

            terminate_processes([proc for proc, _ in self._handles.values()])

            self._handles.clear()
            self._pending_results.clear()

    def _check_running(self) -> None:
        if self._status != _RUNNING:
            raise WorkerError("LocalWorker is not running")

    def __repr__(self) -> str:
        return f"<LocalWorker 0x{id(self):x}>"


class RemoteWorker:
    def __init__(
        self,
        id: str,
        host: str,
        port: int,
        secret: bytes,
        **_kwargs: object,
    ) -> None:
        self.id = id
        self._conn: Any = None
        self._status = _UNINITIALIZED
        self._address = (host, port)
        self._secret = secret
        self._running: dict[int, Node] = {}
        self._threads: int = 0
        self.name = address_to_name((host, port))
        self._log = RemoteAdapter(self.name, logging.getLogger(__name__))

        self._event_handlers: dict[
            tuple[str, str], Callable[[EventType], Iterable[EventType]]
        ] = {
            (_CONNECTING, EVT_HANDSHAKE_RESPONSE): self._event_handshake,
            (_RUNNING, EVT_CAPACITY): self._event_capacity,
            (_RUNNING, EVT_TASK_DONE): self._event_task_done,
            (_RUNNING, EVT_SHUTDOWN): self._event_shutdown,
        }

    @property
    def tasks(self) -> Iterator[Node]:
        yield from self._running.values()

    @property
    def threads(self) -> int:
        return self._threads

    @property
    def handles(self) -> Iterator[HandleType]:
        yield self._conn

    @property
    def status(self) -> str:
        return self._status

    def connect(self, requirements: Collection[Requirement]) -> bool:
        if self._status != _UNINITIALIZED:
            raise WorkerError("Attempted to start already initialized RemoteWorker")

        try:
            self._log.debug("Connecting to worker %s", self.name)
            self._conn = Client(address=self._address, authkey=self._secret)
            self._conn.send(
                {
                    "event": EVT_HANDSHAKE,
                    "cwd": os.getcwd(),
                    "version": paleomix.__version__,
                    "requirements": requirements,
                }
            )
        except (OSError, ProcessError) as error:
            self._log.error("Failed to connect to %s: %s", self.name, error)
            return False
        else:
            self._status = _CONNECTING
            return True

    def start_task(self, task: Node, temp_root: str) -> bool:
        self._check_running()
        self._log.debug("Starting remote task %s with id %s", task, task.id)
        event = {"event": EVT_TASK_START, "task": task, "temp_root": temp_root}
        if not self._send(event):
            return False

        self._running[task.id] = task
        return True

    def has_events(self) -> bool:
        return False

    def get_events(self, _handles: list[HandleType]) -> Iterator[EventType]:
        self._check_running()

        while not self._conn.closed and self._conn.poll():
            try:
                event = self._conn.recv()
            except (ConnectionError, EOFError) as error:
                self._log.error("Connection to worker broke: %s", error)

                if self._status != _TERMINATED:
                    self.shutdown()

                    yield {"event": EVT_SHUTDOWN}
                break

            handler = self._event_handlers.get((self._status, event["event"]))
            if handler is None:
                self._log.error("Unexpected event while %s: %r", self._status, event)
            else:
                self._log.debug("Received %r while %s", event, self._status)

                yield from handler(event)

    def _event_capacity(self, event: EventType) -> Iterable[EventType]:
        self._threads = event["threads"]

        return ()

    def _event_handshake(self, event: EventType) -> Iterator[EventType]:
        if event["error"] is None:
            self._log.debug("Completed handshake")
            self._status = _RUNNING
        else:
            self.shutdown()

        yield event

    def _event_task_done(self, event: EventType) -> Iterator[EventType]:
        event["task"] = self._running.pop(event["task_id"])

        yield event

    def _event_shutdown(self, event: EventType) -> Iterator[EventType]:
        self._status = _TERMINATED
        self._running.clear()
        self._conn.close()

        yield event

    def shutdown(self) -> None:
        self._status = _TERMINATED
        if not self._conn.closed:
            self._log.debug("Shutting down")
            self._send({"event": EVT_SHUTDOWN})
            self._conn.close()

    def _send(self, event: EventType) -> bool:
        self._log.debug("Sending %r", event)
        try:
            self._conn.send(event)
        except OSError as error:
            self._log.error("Failed to send %r: %s", event, error)
            return False
        else:
            return True

    def _check_running(self) -> None:
        if self._status not in (_RUNNING, _CONNECTING):
            raise WorkerError(f"RemoteWorker {self.name} is not running")


def task_wrapper(queue: QueueType, task: Node, temp_root: str) -> None:
    try:
        name = "paleomix task"
        if len(sys.argv) > 1:
            name = f"paleomix {sys.argv[1]} task"

        setproctitle.setproctitle(name)
        # SIGINTs are handled in the main thread only
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        # SIGTERM and SIGHUP are caught to allow proper termination of all sub-commands
        signal.signal(signal.SIGHUP, _task_wrapper_sigterm_handler)
        signal.signal(signal.SIGTERM, _task_wrapper_sigterm_handler)

        task.run(temp_root)
        queue.put((task.id, None, None))
    except NodeError as error:
        backtrace = []
        if error.__cause__ and error.__cause__.__traceback__:
            backtrace = traceback.format_tb(error.__cause__.__traceback__)
            backtrace.append(f"  {error.__cause__!r}")

        queue.put((task.id, error, backtrace))
    except Exception as error:  # noqa: BLE001
        backtrace = traceback.format_tb(error.__traceback__)
        backtrace.append(f"  {error!r}")

        queue.put((task.id, error, backtrace))
    finally:
        queue.close()
        terminate_all_processes()


def _task_wrapper_sigterm_handler(signum: int, _frame: object) -> NoReturn:
    signal.signal(signal.SIGHUP, signal.SIG_DFL)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)

    sys.exit(-signum)
