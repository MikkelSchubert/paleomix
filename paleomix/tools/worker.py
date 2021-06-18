#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import codecs
import json
import logging
import multiprocessing
import os
import random
import signal
import socket
import sys
import uuid
from multiprocessing import Process, ProcessError, Queue, cpu_count
from multiprocessing.connection import Connection, Listener, wait
from typing import Any, Collection, Dict, Iterable, Iterator, List, Optional, Tuple

import paleomix
import paleomix.common.logging
from paleomix.common.argparse import ArgumentParser
from paleomix.common.logging import initialize_console_logging
from paleomix.common.versions import Requirement
from paleomix.core.input import CommandLine, ListTasksEvent, ThreadsEvent
from paleomix.core.workers import (
    AUTO_REGISTER_DIR,
    EVT_CAPACITY,
    EVT_HANDSHAKE,
    EVT_HANDSHAKE_RESPONSE,
    EVT_SHUTDOWN,
    EVT_TASK_DONE,
    EVT_TASK_START,
    QueueType,
    RemoteAdapter,
    _task_wrapper,
    address_to_name,
)
from paleomix.node import Node, NodeError, NodeMissingFilesError
from paleomix.nodegraph import NodeGraph

EventType = Dict[str, Any]
Events = Iterable[EventType]


class Worker:
    def __init__(self, args: argparse.Namespace):
        self._args = args
        self._id: str = args.id
        self._running: Dict[int, Node] = {}
        self._handles: Dict[Any, Tuple[int, Node, Any]] = {}
        self._queue: QueueType = Queue()
        self._address: Tuple[str, int] = (args.host, args.port)
        self._authkey: bytes = args.authkey
        self._threads: int = args.threads
        self._temp_root: str = args.temp_root
        self._filename: Optional[str] = None
        self._interrupted: bool = False
        self._conn: Optional[Connection] = None
        self.name: str = "worker"

        self._commands = {
            EVT_HANDSHAKE: self._manager_handshake,
            EVT_TASK_START: self._manager_start,
        }

        signal.signal(signal.SIGINT, self._sigint_handler)

    def run(self):
        log = logging.getLogger(__name__)
        log.info("Starting worker with %i threads", self._threads)
        if not self._threads:
            log.warning("Worker has no threads; press '+' to allocate more.")

        try:
            listener = Listener(address=self._address, authkey=self._authkey)
        except OSError as error:
            log.error("Could not listen on '%s:%s': %s", *self._address, error)
            return False

        with listener:
            assert isinstance(listener.address, tuple)
            log.info("Listening for PALEOMIX tasks at %s:%s", *listener.address)
            if listener.address[0] in ("127.0.0.1", "127.0.1.1"):
                log.warning("Worker is listening for local connections only!")

            # Create json file containing hostname, port, and secret
            self._filename = self._register(AUTO_REGISTER_DIR, listener.address)

            try:
                self._conn = listener.accept()
            except (OSError, ProcessError) as error:
                log.error("Connection attempt rejected: %s", error)
                return not self._interrupted

            with self._conn:
                address = listener.last_accepted
                assert isinstance(address, tuple)

                name = address_to_name(address)
                self._log = RemoteAdapter(name, logger=log)
                self._log.info("Connection accepted from %s:%i", *address)

                self._main_loop(log)

                return not self._interrupted

    @property
    def tasks(self) -> Iterator[Node]:
        yield from self._running.values()

    def _main_loop(self, log: logging.Logger) -> None:
        assert self._conn is not None

        with CommandLine() as interface:
            while self._running or not self._interrupted:
                handles: List[Any] = [self._conn]
                handles.extend(self._handles)
                handles.extend(interface.handles)

                # Time-out is needed to be able to break on Ctrl+C
                for handle in wait(handles, 1.0):
                    if isinstance(handle, Connection):
                        try:
                            event = handle.recv()
                        except EOFError:
                            log.error("Connection to client broke")
                            return

                        if event["event"] == EVT_SHUTDOWN:
                            self._log.info("Shutting down ..")
                            return

                        func = self._commands.get(event["event"], self._unknown)
                        for event in func(**event):
                            if not self._send(event):
                                return
                    elif handle in interface.handles:
                        if not self._poll_commandline(interface):
                            return
                    elif not self._send(self._join(handle)):
                        return

            self._send({"event": EVT_SHUTDOWN})

    def _register(self, root: str, address: Tuple[str, int]) -> str:
        filename = os.path.join(root, "{}.json".format(uuid.uuid4()))
        host, port = address

        os.makedirs(root, exist_ok=True)
        fd = os.open(filename, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, mode=0o700)
        with open(fd, "w") as handle:
            json.dump(
                {
                    "id": self._id,
                    "host": host,
                    "port": port,
                    "secret": codecs.encode(self._authkey, "base64").decode("utf-8"),
                },
                handle,
            )

        return filename

    def _send(self, event: Dict[str, Any]) -> bool:
        assert self._conn is not None

        self._log.debug("Sending %r", event)
        try:
            self._conn.send(event)
            return True
        except OSError as error:
            self._log.error("Failed to send %r: %s", event, error)
            return False

    def _poll_commandline(self, interface: CommandLine) -> bool:
        for event in interface.process_key_presses():
            if isinstance(event, ListTasksEvent):
                return self._handle_commandline_list_tasks()
            elif isinstance(event, ThreadsEvent):
                return self._handle_commandline_threads(event.change)

        return True

    def _handle_commandline_threads(self, change: int) -> bool:
        threads = min(multiprocessing.cpu_count(), max(0, self._threads + change))
        if threads == self._threads:
            return True

        self._log.info("Max threads changed from %i to %i", self._threads, threads)
        self._threads = threads
        # The number of allocated threads should persist between connections
        self._args.threads = threads

        return self._send({"event": EVT_CAPACITY, "threads": self._threads})

    def _handle_commandline_list_tasks(self):
        tasks = sorted(self.tasks, key=lambda it: it.id)
        threads = sum(task.threads for task in tasks)

        if tasks:
            self._log.info(
                "Running %i tasks on %s (using %i/%i threads):",
                len(tasks),
                self.name,
                threads,
                self._threads,
            )

            for idx, task in enumerate(tasks, start=1):
                self._log.info("  % 2i. %s", idx, task)
        else:
            self._log.info(
                "No tasks running on %s (using 0/%i threads)", self.name, self._threads
            )

        return True

    def _manager_handshake(
        self,
        cwd: str,
        version: str,
        requirements: Collection[Requirement],
        **kwargs: Any,
    ) -> Events:
        if version != paleomix.__version__:
            self._log.error("Version mismatch: %r != %r", version, paleomix.__version__)
            return [
                {
                    "event": EVT_HANDSHAKE_RESPONSE,
                    "error": "Client ({!r}) and worker ({!r}) version mismatch".format(
                        version, paleomix.__version__
                    ),
                }
            ]

        try:
            self._log.info("Changing CWD to %r", cwd)
            os.chdir(cwd)
        except OSError as error:
            self._log.error("Could not change CWD: %r", error)
            return [{"event": EVT_HANDSHAKE_RESPONSE, "error": error}]

        self._log.info("Checking software requirements:")
        if not NodeGraph.check_version_requirements(requirements, force=True):
            return [{"event": EVT_HANDSHAKE_RESPONSE, "error": "Requirements not met"}]

        return [
            {"event": EVT_HANDSHAKE_RESPONSE, "error": None},
            {"event": EVT_CAPACITY, "threads": self._threads},
        ]

    def _manager_start(self, task: Node, temp_root: str, **kwargs: Any) -> Events:
        self._log.info("Starting %s using %i threads", task, task.threads)
        if self._temp_root:
            temp_root = self._temp_root

        proc = Process(
            target=_task_wrapper,
            args=(self._queue, task, temp_root),
            daemon=True,
        )
        proc.start()

        self._handles[proc.sentinel] = (task.id, task, proc)
        self._running[task.id] = task

        threads = sum(task.threads for task in self._running.values())
        if threads > self._threads:
            self._log.warning("Using %i thread(s) too many", threads - self._threads)

        return ()

    def _join(self, handle: Any) -> EventType:
        key, task, proc = self._handles.pop(handle)

        proc.join()
        if proc.exitcode:
            self._log.error("Task process for %s failed with %i", task, proc.exitcode)
            # Task will be the same, but we need to remove it
            task = self._running.pop(key)
            error = NodeError("Process exited with code {}".format(proc.exitcode))
            backtrace = None
        else:
            assert proc.exitcode == 0, proc.exitcode
            # May return result from different task
            key, error, backtrace = self._queue.get()
            task = self._running.pop(key)

            if error is None:
                self._log.info("Finished %s", task)
            elif isinstance(error, NodeMissingFilesError):
                self._log.warning("Finished %s with error %r", task, error)
                self._log.warning(
                    "This may be due to differences in the local/remote filesystem or "
                    "due to file changes not propagating over NFS. If this happens at "
                    "specific filesystem locations, then make sure that the "
                    "local/remote filesystem layout/mount points are identical."
                )
            else:
                self._log.error("Finished %s with error %r", task, error)

        # Signal that the task is done
        return {
            "event": EVT_TASK_DONE,
            "task_id": key,
            "error": error,
            "backtrace": backtrace,
        }

    def _unknown(self, **event: Any):
        self._log.error("Unknown event: %r", event)
        return ()

    def _sigint_handler(self, signum: int, frame: Any) -> None:
        log = logging.getLogger(__name__)
        if self._conn and not self._interrupted:
            self._interrupted = True
            self._threads = 0
            self._send({"event": EVT_CAPACITY, "threads": self._threads})

            log.warning(
                "Keyboard interrupt detected, waiting for running tasks to complete. "
                "Press CTRL-C again to force termination."
            )
        else:
            log.warning("Terminating worker!")
            if self._conn:
                self._conn.close()
            sys.exit(-signum)

    def __enter__(self):
        return self

    def __exit__(self, type: Any, _value: Any, _traceback: Any):
        log = logging.getLogger(__name__)
        for handle, (_, task, proc) in list(self._handles.items()):
            log.warning("Killing %s", task)
            proc.terminate()
            self._join(handle)

        if self._filename and os.path.exists(self._filename):
            log.info("Cleaning up auto-registration JSON at %r", self._filename)
            os.unlink(self._filename)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = ArgumentParser(
        prog="paleomix worker",
        description="Worker process for paleomix pipelines, allowing tasks to be "
        "distributed across multiple systems. Filesystems must be shared between the "
        "pipeline/worker processes using identical layouts/mount-points, including the "
        "home folder which is used to auto-register workers, but software versions are "
        "allowed to differ provided that they meet the basic requirements.",
    )

    group = parser.add_argument_group("Filesystem")
    group.add_argument(
        "--temp-root",
        help="Overrides the temp-root used by the individual pipelines.",
    )

    group = parser.add_argument_group("Worker Scheduling")
    group.add_argument(
        "--host",
        default=socket.gethostname(),
        help="The worker will bind to this host-name",
    )
    group.add_argument(
        "--port",
        type=int,
        default=0,
        help="Port listened on by the worker. Set to 0 to automatically select a port",
    )

    group = parser.add_argument_group("Worker Scheduling")
    group.add_argument(
        "--threads",
        type=int,
        default=cpu_count(),
        help="Maximum number of threads used by worker. Note that the worker may use "
        "more than this number of threads if a task requires a greater number of "
        "threads and no suitable worker is available",
    )
    group.add_argument(
        "--once",
        action="store_true",
        help="Worker will normally serve pipelines indefintely. With this option the "
        "worker will instead terminate after the first pipeline disconnects.",
    )

    paleomix.common.logging.add_argument_group(parser, log_file=False)

    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    args.id = str(uuid.uuid4())

    initialize_console_logging(log_level=args.log_level)

    while True:
        args.authkey = bytes(random.randint(0, 255) for _ in range(64))

        with Worker(args) as worker:
            if not worker.run() or args.once:
                break

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
