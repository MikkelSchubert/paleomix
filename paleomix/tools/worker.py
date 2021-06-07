#!/usr/bin/env python3
# -*- coding: utf8 -*-
import logging
import os
import signal
import socket
import sys

from multiprocessing import Process, Queue, cpu_count
from multiprocessing.connection import Listener, wait, Connection

import paleomix.common.logging

from paleomix.common.argparse import ArgumentParser
from paleomix.common.logging import initialize_console_logging
from paleomix.core.input import CommandLine
from paleomix.core.workers import (
    MANAGER_HANDSHAKE,
    MANAGER_SHUTDOWN,
    MANAGER_START,
    WORKER_FINISHED,
    WORKER_HANDSHAKE,
    WORKER_SHUTDOWN,
    WORKER_THREADS,
    create_secret,
    register_worker,
    _task_wrapper,
    RemoteAdapter,
)
from paleomix.node import NodeError, NodeMissingFilesError
from paleomix.nodegraph import NodeGraph


TERMINATE = object()


class Worker:
    def __init__(self, args):
        self._args = args
        self._running = {}
        self._handles = {}
        self._queue = Queue()
        self._address = (args.host, args.port)
        self._authkey = args.authkey
        self._threads = args.threads
        self._temp_root = args.temp_root
        self._filename = None
        self._interrupted = False
        self._conn = None

        self._commands = {
            MANAGER_HANDSHAKE: self._manager_handshake,
            MANAGER_SHUTDOWN: self._manager_shutdown,
            MANAGER_START: self._manager_start,
        }

        signal.signal(signal.SIGINT, self._sigint_handler)

    def run(self):
        log = logging.getLogger(__name__)
        log.info("Starting worker with %i threads", self._threads)
        if not self._threads:
            log.warning("Worker has no threads; press '+' to allocate more.")

        with Listener(address=self._address, authkey=self._authkey) as listener:
            log.info("Listening for PALEOMIX tasks at %s:%s", *listener.address)
            # Create json file containing hostname, port, and secret
            self._filename = register_worker(*listener.address, self._authkey)

            with listener.accept() as self._conn:
                address = listener.last_accepted
                name = "{}:{}".format(*socket.getnameinfo(address, 0))

                self._log = RemoteAdapter(name, logger=log)
                self._log.info("Connection accepted from %s:%i", *address)

                with CommandLine() as interface:
                    while self._running or not self._interrupted:
                        handles = [self._conn]
                        handles.extend(self._handles)
                        handles.extend(interface.handles)

                        # Time-out is needed to be able to break on Ctrl+C
                        for handle in wait(handles, 5.0):
                            if isinstance(handle, Connection):
                                try:
                                    event = handle.recv()
                                except EOFError:
                                    log.error("Connection to client broke")
                                    return not self._interrupted

                                func = self._commands.get(event["event"], self._unknown)
                                for event in func(**event):
                                    if event is TERMINATE or not self._send(event):
                                        return not self._interrupted
                            elif handle in interface.handles:
                                if not self._poll_commandline(interface):
                                    return not self._interrupted
                            elif not self._send(self._join(handle)):
                                return not self._interrupted

                    # Report the termination of any remaining running tasks
                    self._send({"event": WORKER_SHUTDOWN})

        return not self._interrupted

    def _send(self, event):
        self._log.debug("Sending %r", event)
        try:
            self._conn.send(event)
            return True
        except OSError as error:
            self._log.error("Failed to send %r: %s", event, error)
            return False

    def _poll_commandline(self, interface):
        new_threads = interface.process_key_presses(
            threads=self._threads,
            tasks=self._running.values(),
        )

        if self._threads == new_threads:
            return True

        self._threads = new_threads
        # The number of allocated threads should persist between connections
        self._args.threads = new_threads

        return self._send({"event": WORKER_THREADS, "threads": self._threads})

    def _manager_handshake(self, cwd, requirements, **kwargs):
        try:
            self._log.info("Changing CWD to %r", cwd)
            os.chdir(cwd)
        except OSError as error:
            self._log.error("Could not change CWD: %r", error)
            return [{"event": WORKER_HANDSHAKE, "error": error}]

        self._log.info("Checking software requirements:")
        if not NodeGraph.check_version_requirements(requirements, force=True):
            return [{"event": WORKER_HANDSHAKE, "error": "Requirements not met"}]

        return [
            {"event": WORKER_HANDSHAKE, "error": None},
            {"event": WORKER_THREADS, "threads": self._threads},
        ]

    def _manager_shutdown(self, **kwargs):
        self._log.info("Shutting down ..")

        return [TERMINATE]

    def _manager_start(self, key, task, temp_root, **kwargs):
        self._log.info("Starting %s using %i threads", task, task.threads)
        if self._temp_root:
            temp_root = self._temp_root

        proc = Process(
            target=_task_wrapper,
            args=(self._queue, key, task, temp_root),
            daemon=True,
        )
        proc.start()

        self._handles[proc.sentinel] = (key, task, proc)
        self._running[key] = task

        threads = sum(task.threads for task in self._running.values())
        if threads > self._threads:
            self._log.warning("Using %i thread(s) too many", threads - self._threads)

        return ()

    def _join(self, handle):
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
                self._log.info("Finished task %s", task)
            elif isinstance(error, NodeMissingFilesError):
                self._log.warning("Finished task %s with error %r", task, error)
                self._log.warning(
                    "This may be due to differences in the local/remote filesystem or "
                    "due to file changes not propagating over NFS. If this happens at "
                    "specific filesystem locations, then make sure that the "
                    "local/remote filesystem layout/mount points are identical."
                )
            else:
                self._log.error("Finished task %s with error %r", task, error)

        # Signal that the task is done
        return {
            "event": WORKER_FINISHED,
            "key": key,
            "error": error,
            "backtrace": backtrace,
        }

    def _unknown(self, **event):
        self._log.error("Unknown event: %r", event)
        return ()

    def _sigint_handler(self, signum, frame):
        log = logging.getLogger(__name__)
        if self._conn and not self._interrupted:
            self._interrupted = True
            self._threads = 0
            self._send({"event": WORKER_THREADS, "threads": self._threads})

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

    def __exit__(self, type, _value, _traceback):
        log = logging.getLogger(__name__)
        for handle, (_key, task, proc) in list(self._handles.items()):
            log.warning("Killing %s", task)
            proc.terminate()
            self._join(handle)

        if self._filename and os.path.exists(self._filename):
            log.info("Cleaning up auto-registration JSON at %r", self._filename)
            os.unlink(self._filename)


def parse_args(argv):
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


def main(argv):
    args = parse_args(argv)
    initialize_console_logging(log_level=args.log_level)

    while True:
        args.authkey = create_secret()

        with Worker(args) as worker:
            if not worker.run() or args.once:
                break

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
