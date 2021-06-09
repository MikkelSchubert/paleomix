import codecs
import json
import logging
import multiprocessing
import os
import os.path
import random
import signal
import socket
import sys
import time
import traceback
import uuid

from multiprocessing.connection import Client, wait

import paleomix.common.system

from paleomix.core.input import CommandLine
from paleomix.nodegraph import NodeGraph
from paleomix.node import NodeError


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


#
_UNINITIALIZED = "uninitialized"
_CONNECTING = "connecting"
_RUNNING = "running"
_TERMINATED = "terminated"


class WorkerError(RuntimeError):
    pass


class Manager:
    def __init__(self, threads, requirements, temp_root, auto_connect=True):
        self._local = None
        self._interface = CommandLine()
        self._workers = {}
        self._requirements = requirements
        self._threads = threads
        self._temp_root = temp_root
        self._next_auto_connect = 0 if auto_connect else float("inf")
        self._log = logging.getLogger(__name__)

    @property
    def workers(self):
        return dict(self._workers)

    @property
    def tasks(self):
        for worker in self._workers.values():
            yield from worker.tasks

    def start(self):
        if self._local is not None:
            raise WorkerError("Manager already started")

        local_worker = LocalWorker(self, self._interface, self._threads)
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

    def wait_for_workers(self):
        self._check_started()
        while not all(worker.status != _RUNNING for worker in self._workers.values()):
            for event in self.poll():
                if event["event"] == EVT_HANDSHAKE_RESPONSE:
                    if event["error"] is not None:
                        return False

        return True

    def shutdown(self):
        for worker in self._workers.values():
            worker.shutdown()
        self._workers.clear()
        self._local = None

    def poll(self, blocking=False):
        self._check_started()
        self._auto_connect_to_workers(blocking)

        for worker, event in self._collect_events():
            if not (isinstance(event, dict) and "event" in event):
                self._log.error("Malformed event from %s: %r", worker.name, event)
                continue

            event["worker"] = worker.id
            event["worker_name"] = worker.name

            if event["event"] == EVT_HANDSHAKE_RESPONSE:
                if event["error"] is not None:
                    self._log.error("Handshake with worker %r failed", worker.name)
                    self._workers.pop(worker.id)

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

    def start_task(self, worker_id, task):
        self._check_started()

        worker = self._workers.get(worker_id)
        if worker is None:
            self._log.error("Tried to start task on unknown worker %r", worker_id)
            return False

        return worker.start_task(task, self._temp_root)

    def _collect_events(self):
        events = []
        timeout = 5.0
        while self._local and self._local.events:
            for event in self._local.events:
                events.append((self._local, event))

            self._local.events.clear()
            # No need to block if we already have events to act on
            timeout = 0

        handles = {}
        for worker in self._workers.values():
            for handle in worker.handles:
                handles[handle] = worker

        for handle in wait(handles, timeout):
            worker = handles[handle]
            for event in worker.get(handle):
                events.append((worker, event))

        return events

    def _auto_connect_to_workers(self, interval=15.0):
        any_errors = False
        if self._next_auto_connect <= time.monotonic():
            for info in collect_workers():
                self._log.info("Connecting to %s:%s", info["host"], info["port"])
                worker = RemoteRunner(**info)
                if worker.id in self._workers:
                    self._log.error("Already connected to worker with id", worker.id)
                    any_errors = True
                elif worker.connect(self._requirements):
                    self._workers[worker.id] = worker
                else:
                    any_errors = True

            self._next_auto_connect = time.monotonic() + interval

        return not any_errors

    def _check_started(self):
        if self._local is None:
            raise WorkerError("Manager not started")

    def __enter__(self):
        self._interface.setup()
        return self

    def __exit__(self, _type, _value, _traceback):
        self._interface.teardown()
        self.shutdown()


class RemoteAdapter(logging.LoggerAdapter):
    def __init__(self, name, logger):
        super().__init__(logger, {"remote": name})

    def process(self, msg, kwargs):
        return "[{}] {}".format(self.extra["remote"], msg), kwargs


class LocalWorker:
    def __init__(self, manager, interface, threads):
        self.id = str(uuid.uuid4())
        self._manager = manager
        self._interface = interface
        self._threads = threads
        self._queue = multiprocessing.Queue()
        self._handles = {}
        self._running = {}
        self._status = _UNINITIALIZED
        self._log = logging.getLogger(__name__)
        self.name = "localhost"
        self.events = []

    @property
    def tasks(self):
        yield from self._running.values()

    @property
    def threads(self):
        return self._threads

    @property
    def handles(self):
        handles = list(self._handles)
        handles.extend(self._interface.handles)

        return handles

    @property
    def status(self):
        return self._status

    def connect(self, requirements):
        if self._status != _UNINITIALIZED:
            raise WorkerError("Attempted to start already initialized LocalWorker")

        self._log.info("Checking required software on localhost")
        if not NodeGraph.check_version_requirements(requirements):
            return False

        self._status = _RUNNING
        self.events.append({"event": EVT_HANDSHAKE_RESPONSE, "error": None})

        return True

    def start_task(self, task, temp_root):
        self._check_running()

        self._log.debug("Starting local task %s with id %s", task, task.id)
        proc = multiprocessing.Process(
            target=_task_wrapper,
            args=(self._queue, task, temp_root),
            daemon=True,
        )

        proc.start()
        self._handles[proc.sentinel] = (proc, task)
        self._running[task.id] = task

        return True

    def get(self, handle):
        self._check_running()

        if handle in self._interface.handles:
            self._threads = self._interface.process_key_presses(
                threads=self._threads,
                workers=self._manager.workers.values(),
            )

            return ()

        # The process for this task has completed, but the results may already have been
        # returned to the pipeline in an earlier call, since it is returned via a queue.
        proc, task = self._handles.pop(handle)

        proc.join()
        if proc.exitcode:
            self._log.debug("Local join of task %s with failed", task)
            message = "Process terminated with exit code {}".format(proc.exitcode)

            error = NodeError(message)
            backtrace = None
        else:
            self._log.debug("Joined local task %s", task)

            # This can return results for a different task than the one joined above
            key, error, backtrace = self._queue.get()
            task = self._running.pop(key)

        # Signal that the task is done
        return [
            {
                "event": EVT_TASK_DONE,
                "task": task,
                "error": error,
                "backtrace": backtrace,
            },
        ]

    def shutdown(self):
        if self._status != _TERMINATED:
            self._status = _TERMINATED
            self._log.debug("Shutting down local worker")
            for proc, task in tuple(self._handles.values()):
                self._log.warning("Killing task %s", task)
                proc.terminate()
            self._handles.clear()
            self._running.clear()

    def _check_running(self):
        if self._status != _RUNNING:
            raise WorkerError("LocalWorker is not running")


class RemoteRunner:
    def __init__(self, id, host, port, secret, **_kwargs):
        self.id = id
        self._conn = None
        self._status = _UNINITIALIZED
        self._address = (host, port)
        self._secret = secret
        self._running = {}
        self._threads = 0
        self.name = "{}:{}".format(*socket.getnameinfo((host, port), 0))
        self._log = RemoteAdapter(self.name, logging.getLogger(__name__))

        self._event_handlers = {
            (_CONNECTING, EVT_HANDSHAKE_RESPONSE): self._event_handshake,
            (_RUNNING, EVT_CAPACITY): self._event_capacity,
            (_RUNNING, EVT_TASK_DONE): self._event_task_done,
            (_RUNNING, EVT_SHUTDOWN): self._event_shutdown,
        }

    @property
    def tasks(self):
        yield from self._running.values()

    @property
    def threads(self):
        return self._threads

    @property
    def handles(self):
        yield self._conn

    @property
    def status(self):
        return self._status

    def connect(self, requirements):
        if self._status != _UNINITIALIZED:
            raise WorkerError("Attempted to start already initialized RemoteWorker")

        try:
            self._log.debug("Connecting to worker %s", self.name)
            self._conn = Client(address=self._address, authkey=self._secret)
            self._conn.send(
                {
                    "event": EVT_HANDSHAKE,
                    "cwd": os.getcwd(),
                    "requirements": requirements,
                }
            )

            self._status = _CONNECTING

            return True
        except OSError as error:
            self._log.error("Failed to connect to %s: %s", self.name, error)
            return False

    def start_task(self, task, temp_root):
        self._check_running()
        self._log.debug("Starting remote task %s with id %s", task, task.id)
        event = {"event": EVT_TASK_START, "task": task, "temp_root": temp_root}
        if not self._send(event):
            return False

        self._running[task.id] = task
        return True

    def get(self, handle):
        self._check_running()
        if handle is not self._conn:
            raise ValueError(handle)

        while not self._conn.closed and self._conn.poll():
            try:
                event = self._conn.recv()
            except EOFError:
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

    def _event_capacity(self, event):
        self._threads = event["threads"]

        return ()

    def _event_handshake(self, event):
        if event["error"] is None:
            self._log.debug("Completed handshake")
            self._status = _RUNNING

        yield event

    def _event_task_done(self, event):
        event["task"] = self._running.pop(event["task_id"])

        yield event

    def _event_shutdown(self, event):
        self._status = _TERMINATED
        self._running.clear()
        self._conn.close()

        yield event

    def shutdown(self):
        self._status = _TERMINATED
        if not self._conn.closed:
            self._log.debug("Shutting down")
            self._send({"event": EVT_SHUTDOWN})
            self._conn.close()

    def _send(self, event):
        self._log.debug("Sending %r", event)
        try:
            self._conn.send(event)
            return True
        except OSError as error:
            self._log.error("Failed to send %r: %s", event, error)
            return False

    def _check_running(self):
        if self._status not in (_RUNNING, _CONNECTING):
            raise WorkerError("RemoteWorker {} is not running".format(self.name))


def register_worker(id, host, port, secret, root="~/.paleomix/remote/"):
    root = os.path.expanduser(root)
    os.makedirs(root, exist_ok=True)
    filename = os.path.join(root, "{}.json".format(uuid.uuid4()))

    fd = os.open(filename, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, mode=0o700)
    with open(fd, "w") as handle:
        json.dump(
            {
                "id": id,
                "host": host,
                "port": port,
                "secret": codecs.encode(secret, "base64").decode("utf-8"),
            },
            handle,
        )

    return filename


def collect_workers(root="~/.paleomix/remote/"):
    root = os.path.expanduser(root)
    for filename in os.listdir(root):
        filename = os.path.join(root, filename)

        if filename.endswith(".json") and os.path.isfile(filename):
            with open(filename, "rt") as handle:
                try:
                    data = json.load(handle)
                except json.decoder.JSONDecodeError:
                    continue
            os.unlink(filename)

            data["secret"] = codecs.decode(data["secret"].encode("utf-8"), "base64")

            yield data


def create_secret():
    return bytes(random.randint(0, 255) for _ in range(64))


def _task_wrapper(queue, task, temp_root):
    name = "paleomix task"
    if len(sys.argv) > 1:
        name = "paleomix {} task".format(sys.argv[1])

    paleomix.common.system.set_procname(name)
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    try:
        task.run(temp_root)
        queue.put((task.id, None, None))
    except NodeError as error:
        backtrace = []
        if error.__cause__ and error.__cause__.__traceback__:
            backtrace = traceback.format_tb(error.__cause__.__traceback__)
            backtrace.append("  {!r}".format(error.__cause__))

        queue.put((task.id, error, backtrace))
    except Exception:
        _exc_type, exc_value, exc_traceback = sys.exc_info()
        backtrace = traceback.format_tb(exc_traceback)
        backtrace.append("  {!r}".format(exc_value))

        queue.put((task.id, exc_value, backtrace))
