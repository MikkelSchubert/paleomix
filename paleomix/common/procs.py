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
from __future__ import annotations

import atexit
import os
import shlex
import signal
import sys
from os import PathLike
from subprocess import Popen, TimeoutExpired
from typing import IO, Any, AnyStr, Iterable, List, Optional, cast

# List of running processes; can be terminated with `terminate_all_processes`
_RUNNING_PROCS: List[Popen[Any]] = []


def quote_args(args: Any) -> str:
    if isinstance(args, (str, bytes, os.PathLike)):
        args = [args]

    values: List[str] = []
    for arg in args:
        if isinstance(arg, os.PathLike):
            arg = os.fspath(cast(os.PathLike[AnyStr], arg))

        if isinstance(arg, bytes):
            arg = arg.decode("utf-8", errors="replace")

        values.append(shlex.quote(arg))

    return " ".join(values)


def join_procs(procs: Iterable[Popen[Any]], out: IO[str] = sys.stderr):
    """Joins a set of Popen processes. If a processes fail, the remaining
    processes are terminated. The function returns a list of return-code,
    containing the result of each call. Status messages are written to STDERR
    by default.
    """
    sleep_time = 0.05
    commands = list(enumerate(procs))
    return_codes: List[Optional[int]] = [None] * len(commands)

    assert all(hasattr(cmd, "args") for (_, cmd) in commands)

    print("Joining subprocesses:", file=out)
    while commands and not any(return_codes):
        try:
            # Wait for arbitrary command
            commands[0][1].wait(sleep_time if len(commands) > 1 else None)
        except TimeoutExpired:
            sleep_time = min(1, sleep_time * 2)

        for (index, command) in list(commands):
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


def register_process(proc: Popen[Any]) -> None:
    """Register a process for automatic/forced termination."""
    _RUNNING_PROCS.append(proc)


def unregister_process(proc: Popen[Any]) -> None:
    """Unregister a process for automatic/forced termination."""
    if proc in _RUNNING_PROCS:
        _RUNNING_PROCS.remove(proc)


@atexit.register
def terminate_all_processes() -> None:
    """Terminate all registered processes. Must be called in signal handlers."""
    while _RUNNING_PROCS:
        proc = _RUNNING_PROCS.pop()

        try:
            proc.terminate()
        except OSError:
            # Ignore already closed processes, etc.
            pass
