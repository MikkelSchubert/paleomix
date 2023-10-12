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

import atexit
import contextlib
import os
import shlex
import signal
import sys
from subprocess import Popen, TimeoutExpired
from typing import IO, TYPE_CHECKING, Any, Iterable, cast

if TYPE_CHECKING:
    from typing_extensions import Protocol

    class _SupportsTerminate(Protocol):
        def terminate(self) -> None:
            ...


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
_RUNNING_PROCS: list[_SupportsTerminate] = []


def running_processes() -> list[_SupportsTerminate]:
    return list(_RUNNING_PROCS)


def register_process(proc: _SupportsTerminate) -> None:
    """Register a process for automatic/forced termination."""
    _RUNNING_PROCS.append(proc)


def unregister_process(proc: _SupportsTerminate) -> None:
    """Unregister a process for automatic/forced termination."""
    if proc in _RUNNING_PROCS:
        _RUNNING_PROCS.remove(proc)


@atexit.register
def terminate_all_processes() -> None:
    """Terminate all registered processes. Must be called in signal handlers."""
    while _RUNNING_PROCS:
        proc = _RUNNING_PROCS.pop()

        with contextlib.suppress(OSError):
            # Ignore already closed processes, etc.
            proc.terminate()
