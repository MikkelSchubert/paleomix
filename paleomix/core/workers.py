import codecs
import json
import logging
import multiprocessing
import os
import os.path
import random
import signal
import sys
import time
import traceback
import uuid

from multiprocessing.connection import Client, wait

import paleomix.common.system

from paleomix.core.input import CommandLine
from paleomix.nodegraph import NodeGraph
from paleomix.node import NodeError


# Messages sent by the manager (pipeline):
#
MANAGER_HANDSHAKE = "M:HANDSHAKE"
MANAGER_SHUTDOWN = "M:SHUTDOWN"
MANAGER_START = "M:START"

# Messages sent by the worker
WORKER_HANDSHAKE = "W:HANDSHAKE"
WORKER_THREADS = "W:THREADS"
WORKER_SHUTDOWN = "W:SHUTDOWN"
WORKER_FINISHED = "W:FINISHED"


UNINITIALIZED = 0
RUNNING = 1
TERMINATED = 2


class WorkerError(RuntimeError):
    pass


class Manager:
    def __init__(self, threads, requirements, temp_root, auto_connect=True):
        self._local = None
        self._interface = CommandLine()
        self._workers = []
        self._requirements = requirements
        self._threads = threads
        self._temp_root = temp_root
        self._auto_connect = auto_connect
        self._next_auto_connect = 0
        self._log = logging.getLogger(__name__)

    def connect(self):
        if self._workers or self._local:
            raise WorkerError("Manager already started")

        self._local = LocalWorker(self, self._interface, self._threads)
        if not self._local.connect(self._requirements):
            # The local worker failing is a fatal error since it cannot be fixed
            # without restarting the entire pipeline
            return False

        self._workers.append(self._local)

        if self._threads <= 0:
            self._log.warning(
                "Local worker process has no threads assigned; either "
                "increase allocation with '+' or start worker processes."
            )

        return self._auto_connect_to_workers()

    def shutdown(self):
        for worker in self._workers:
            worker.shutdown()

    def poll(self):
        self._auto_connect_to_workers(blocking=not self._workers)

        handles = {}
        for worker in self._workers:
            for handle in worker.handles:
                handles[handle] = worker

        for handle in wait(handles, 0.0 if self.idle else 5.0):
            worker = handles[handle]

            for event in worker.get(handle):
                event["worker"] = worker
                # Worker terminated gracefully, but may have had still running tasks
                if event["event"] == WORKER_SHUTDOWN:
                    self._workers.remove(worker)
                    worker.shutdown()

                yield event

        max_threads = max((worker.threads for worker in self._workers), default=0)
        for worker in self._workers:
            idle_threads = worker.idle_threads
            if idle_threads > 0:
                # Signal that there are available tasks
                yield {
                    "event": WORKER_THREADS,
                    "threads": idle_threads,
                    # The biggest workers are allowed to exceed capacity if needed
                    "overcommit": max_threads and idle_threads == max_threads,
                    "worker": worker,
                }

    def start(self, worker, task):
        if not worker.start(task, self._temp_root):
            return False

        return True

    @property
    def ready(self):
        for worker in self._workers:
            if worker.status != RUNNING:
                return False

        return True

    @property
    def idle(self):
        return not any(self.tasks)

    @property
    def tasks(self):
        for worker in self._workers:
            yield from worker.tasks

    def _auto_connect_to_workers(self, blocking=False, interval=15.0):
        any_errors = False
        while self._auto_connect:
            if self._next_auto_connect <= time.monotonic():
                for info in collect_workers():
                    self._log.info("Connecting to %s:%s", info["host"], info["port"])
                    worker = RemoteRunner(info["host"], info["port"], info["secret"])
                    if not worker.connect(self._requirements):
                        any_errors = True

                    self._workers.append(worker)

                self._next_auto_connect = time.monotonic() + 15

            if self._workers or not blocking:
                break

            self._log.debug("Waiting for workers to be registered ..")
            time.sleep(self._next_auto_connect - time.monotonic())

        return not any_errors

    def __enter__(self):
        self._interface.setup()
        return self

    def __exit__(self, _type, _value, _traceback):
        self._interface.teardown()


class RemoteAdapter(logging.LoggerAdapter):
    def __init__(self, name, logger):
        super().__init__(logger, {"remote": name})

    def process(self, msg, kwargs):
        return "[{}] {}".format(self.extra["remote"], msg), kwargs


class LocalWorker:
    def __init__(self, manager, interface, threads):
        self._manager = manager
        self._interface = interface
        self._threads = threads
        self._queue = multiprocessing.Queue()
        self._handles = {}
        self._running = {}
        self._status = UNINITIALIZED
        self._log = logging.getLogger(__name__)
        self.name = ":builtin:"

    @property
    def tasks(self):
        yield from self._running.values()

    @property
    def threads(self):
        return self._threads

    @property
    def idle_threads(self):
        return self._threads - sum(task.threads for task in self._running.values())

    @property
    def handles(self):
        handles = list(self._handles.keys())
        handles.extend(self._interface.handles)

        return handles

    @property
    def status(self):
        return self._status

    def connect(self, requirements):
        self._log.info("Checking required software on localhost")
        if not NodeGraph.check_version_requirements(requirements):
            self._log.error(
                "Please refer to the PALEOMIX installation instructions at\n"
                "  https://paleomix.readthedocs.io/en/stable/"
            )

        self._status = RUNNING
        return True

    def start(self, task, temp_root):
        key = str(uuid.uuid4())
        proc = multiprocessing.Process(
            target=_task_wrapper,
            args=(self._queue, key, task, temp_root),
            daemon=True,
        )

        proc.start()
        self._log.debug("Started local task %s with key %s", task, key)

        self._handles[proc.sentinel] = (proc, task)
        self._running[key] = task

        return True

    def get(self, handle):
        if handle in self._interface.handles:
            self._threads = self._interface.process_key_presses(
                threads=self._threads,
                tasks=self._manager.tasks,
            )

            return ()

        proc, task = self._handles.pop(handle)

        proc.join()
        if proc.exitcode:
            self._log.debug("Local join of %s with failed", task)
            message = "Process terminated with exit code {}".format(proc.exitcode)

            error = NodeError(message)
            backtrace = None
        else:
            assert proc.exitcode == 0, proc.exitcode
            self._log.debug("Joined local task %s", task)

            # This can return results for a different task than the one joined
            key, error, backtrace = self._queue.get()
            task = self._running.pop(key)

        # Signal that the task is done
        return [
            {
                "event": WORKER_FINISHED,
                "task": task,
                "error": error,
                "backtrace": backtrace,
            }
        ]

    def shutdown(self):
        self._status = TERMINATED
        self._log.debug("Shutting down local worker")
        for proc, task in tuple(self._handles.values()):
            self._log.warning("Killing task %s", task)
            proc.terminate()


class RemoteRunner:
    def __init__(self, host, port, secret):
        self._conn = None
        self._status = UNINITIALIZED
        self._address = (host, port)
        self._secret = secret
        self._running = {}
        self._threads = 0
        self.name = "{}:{}".format(host, port)
        self._log = RemoteAdapter(self.name, logging.getLogger(__name__))

    @property
    def tasks(self):
        yield from self._running.values()

    @property
    def threads(self):
        return self._threads

    @property
    def idle_threads(self):
        return self._threads - sum(task.threads for task in self._running.values())

    @property
    def handles(self):
        yield self._conn

    @property
    def status(self):
        return self._status

    def connect(self, requirements):
        self._log.debug("Connecting to worker %s", self.name)

        try:
            self._conn = Client(address=self._address, authkey=self._secret)
            self._conn.send(
                {
                    "event": MANAGER_HANDSHAKE,
                    "cwd": os.getcwd(),
                    "requirements": requirements,
                }
            )

            return True
        except OSError as error:
            self._log.error("Failed to connect to %s: %s", self.name, error)
            return False

    def start(self, task, temp_root):
        if self._status != RUNNING:
            raise RuntimeError("%r used while not ready" % (self._address))

        key = str(uuid.uuid4())
        self._log.debug("Starting remote task %s with key %s", task, key)

        event = {
            "event": MANAGER_START,
            "key": key,
            "task": task,
            "temp_root": temp_root,
        }

        if not self._send(event):
            return False

        self._running[key] = task
        return True

    def get(self, handle):
        if handle is not self._conn:
            raise ValueError(handle)

        while not self._conn.closed and self._conn.poll():
            try:
                event = self._conn.recv()
            except EOFError:
                if self._status != TERMINATED:
                    self._status = TERMINATED
                    tasks = list(self._running.values())
                    self._running.clear()
                    yield {"event": WORKER_SHUTDOWN, "tasks": tasks}
                break

            self._log.debug("Received %r", event)
            if self._status == RUNNING:
                if event["event"] == WORKER_THREADS:
                    self._threads = event["threads"]
                elif event["event"] == WORKER_FINISHED:
                    # Signal that the task is done
                    yield {
                        "event": WORKER_FINISHED,
                        "task": self._running.pop(event["key"]),
                        "error": event["error"],
                        "backtrace": event["backtrace"],
                    }
                elif event["event"] == WORKER_SHUTDOWN:
                    self._status = TERMINATED
                    tasks = list(self._running.values())
                    self._running.clear()
                    yield {"event": WORKER_SHUTDOWN, "tasks": tasks}
                else:
                    self._log.error("Unknown event %r", event)
            elif self._status == UNINITIALIZED:
                if event["event"] == WORKER_HANDSHAKE:
                    if event["error"] is None:
                        self._log.debug("Completed handshake")
                        self._status = RUNNING
                    else:
                        self._log.error("Handshake failed: %s", event["error"])
                else:
                    self._log.error("Got message while NOT ready: %r", event)
            else:
                self._log.error("Got message while terminated: %r", event)

    def shutdown(self):
        if not self._conn.closed:
            self._log.debug("Shutting down")
            self._send({"event": MANAGER_SHUTDOWN})
            self._conn.close()

    def _send(self, event):
        self._log.debug("Sending %r", event)
        try:
            self._conn.send(event)
            return True
        except OSError as error:
            self._log.error("Failed to send %r: %s", event, error)
            return False


def register_worker(host, port, secret, root="~/.paleomix/remote/"):
    root = os.path.expanduser(root)
    os.makedirs(root, exist_ok=True)
    filename = os.path.join(root, "{}.json".format(uuid.uuid4()))

    fd = os.open(filename, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, mode=0o700)
    with open(fd, "w") as handle:
        json.dump(
            {
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


def _task_wrapper(queue, key, task, temp_root):
    name = "paleomix task"
    if len(sys.argv) > 1:
        name = "paleomix {} task".format(sys.argv[1])

    paleomix.common.system.set_procname(name)
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    try:
        task.run(temp_root)
        queue.put((key, None, None))
    except NodeError as error:
        backtrace = []
        if error.__cause__ and error.__cause__.__traceback__:
            backtrace = traceback.format_tb(error.__cause__.__traceback__)
            backtrace.append("  {!r}".format(error.__cause__))

        queue.put((key, error, backtrace))
    except Exception:
        _exc_type, exc_value, exc_traceback = sys.exc_info()
        backtrace = traceback.format_tb(exc_traceback)
        backtrace.append("  {!r}".format(exc_value))

        queue.put((key, exc_value, backtrace))
