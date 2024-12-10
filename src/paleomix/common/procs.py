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

import contextlib
import logging
import os
import shlex
import signal
import subprocess
import sys
import time
from collections import defaultdict
from multiprocessing import Process
from pathlib import Path
from subprocess import Popen, TimeoutExpired
from typing import IO, TYPE_CHECKING, Any, Iterable, Sequence, Union, cast

from typing_extensions import TypeAlias

ProcessTypes: TypeAlias = Union["Popen[str]", "Popen[bytes]", Process]


PIPE = subprocess.PIPE


def quote_args(args: object) -> str:
    objects: Iterable[object]
    if isinstance(args, os.PathLike):
        # WORKAROUND for pyright warning about "partially unknown" types
        objects = [cast(object, args)]
    elif isinstance(args, (str, bytes)) or not isinstance(args, Iterable):
        objects = [args]
    else:
        objects = args

    values: list[str] = []
    for value in objects:
        if isinstance(value, os.PathLike):
            value = os.fsdecode(cast(Any, value))

        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")

        values.append(shlex.quote(str(value)))

    return " ".join(values)


def join_procs(
    procs: Iterable[Popen[Any]],
    out: IO[str] = sys.stderr,
) -> list[int | None]:
    """Joins a set of Popen processes. If a processes fail, the remaining
    processes are terminated. The function returns a list of return-code,
    containing the result of each call. Status messages are written to STDERR
    by default.
    """
    sleep_time = 0.05
    commands = list(enumerate(procs))
    return_codes: list[int | None] = [None] * len(commands)

    print("Joining subprocesses:", file=out)
    while commands and not any(return_codes):
        try:
            # Wait for arbitrary command
            commands[0][1].wait(sleep_time if len(commands) > 1 else None)
        except TimeoutExpired:
            sleep_time = min(1, sleep_time * 2)

        for index, command in list(commands):
            if command.poll() is not None:
                return_code = command.wait()
                return_codes[index] = return_code
                commands.remove((index, command))
                sleep_time = 0.05

                if return_code < 0:
                    return_code = signal.Signals(-return_code).name

                print(f"  - Command finished: {quote_args(command.args)}", file=out)
                print(f"    Return-code:      {return_code}", file=out)

    if any(return_codes):
        for index, command in commands:
            print(f"  - Terminating command: {quote_args(command.args)}", file=out)
            command.terminate()
            return_codes[index] = command.wait()

        print("Errors occurred during processing!", file=out)

    return return_codes


# List of running processes; can be terminated with `terminate_all_processes`
_RUNNING_PROCS: dict[int, list[ProcessTypes]] = defaultdict(list)


class RegisteredProcess(Process):
    def start(self) -> None:
        _register_process(self)
        super().start()

    def join(self, timeout: float | None = None) -> None:
        super().join(timeout)
        _unregister_process(self)


PopenBase = Popen[bytes] if TYPE_CHECKING else Popen


class RegisteredPopen(PopenBase):
    def __init__(
        self,
        args: str | bytes | Sequence[str | bytes],
        *,
        stdin: None | int | IO[bytes] = None,
        stdout: None | int | IO[bytes] = None,
        stderr: None | int | IO[bytes] = None,
        close_fds: bool = True,
        cwd: str | Path | None = None,
        start_new_session: bool = False,
    ) -> None:
        super().__init__(
            args=args,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            close_fds=close_fds,
            cwd=cwd,
            start_new_session=start_new_session,
        )

        _register_process(self)

    def wait(self, timeout: float | None = None) -> int:
        return_code = super().wait(timeout)
        _unregister_process(self)

        return return_code


def running_processes() -> list[ProcessTypes]:
    return list(_RUNNING_PROCS[os.getpid()])


def terminate_processes(
    processes: Iterable[ProcessTypes],
    timeout: float | None = None,
) -> None:
    log = logging.getLogger(__name__)
    start = time.time()

    processes = tuple(processes)
    for proc in processes:
        # Ignore already closed processes, etc.
        with contextlib.suppress(OSError):
            if isinstance(proc, Popen):
                command = quote_args(proc.args)
                if len(command) > 80:
                    command = command[:77] + "..."

                log.warning("Terminating process %s: %s", proc.pid, command)
            else:
                log.warning("Terminating worker processes %s", proc.pid)

            proc.terminate()

    for proc in processes:
        time_left = None
        if timeout is not None:
            time_left = max(0.0, timeout - (time.time() - start))

        with contextlib.suppress(TimeoutError):
            if isinstance(proc, Process):
                proc.join(timeout=time_left)
            else:
                proc.wait(timeout=time_left)

            log.debug("Joined %s (timeout=%r)", proc, time_left)


def terminate_all_processes(timeout: float | None = None) -> None:
    """Terminate all registered processes. Must be called in signal handlers."""
    terminate_processes(running_processes(), timeout=timeout)


def _register_process(proc: ProcessTypes) -> None:
    """Register a process for automatic/forced termination."""
    log = logging.getLogger(__name__)
    log.debug("[%s] register :%s", os.getpid(), proc)

    _RUNNING_PROCS[os.getpid()].append(proc)


def _unregister_process(proc: ProcessTypes) -> None:
    """Unregister a process for automatic/forced termination."""
    processes = _RUNNING_PROCS[os.getpid()]
    # May have been removed by `terminate_all_processes`, exception handlers, etc.
    if proc in processes:
        log = logging.getLogger(__name__)
        log.debug("[%s] unregister %s", os.getpid(), proc)

        processes.remove(proc)
