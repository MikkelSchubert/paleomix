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

import signal
import sys
import time
from subprocess import Popen
from typing import IO, Any, Iterable, List, Optional, cast


def join_procs(procs: Iterable[Popen[Any]], out: IO[str] = sys.stderr):
    """Joins a set of Popen processes. If a processes fail, the remaining
    processes are terminated. The function returns a list of return-code,
    containing the result of each call. Status messages are written to STDERR
    by default.
    """
    sleep_time = 0.05
    commands = list(enumerate(procs))
    assert all(hasattr(cmd, "args") for (_, cmd) in commands)

    return_codes = [None] * len(commands)  # type: List[Optional[int]]
    out.write("Joinining subprocesses:\n")
    while commands:
        for (index, command) in list(commands):
            if command.poll() is not None:
                return_code = command.wait()
                return_codes[index] = return_code
                commands.remove((index, command))
                sleep_time = 0.05

                if return_code < 0:
                    return_code = signal.Signals(-return_code).name

                out.write(
                    "  - Command finished: %s\n"
                    "    Return-code:      %s\n"
                    % (" ".join(cast(Any, command).args), return_code)
                )
                out.flush()
            elif any(return_codes):
                out.write(
                    "  - Terminating command: %s\n"
                    % (" ".join(cast(Any, command).args),)
                )
                out.flush()

                command.terminate()
                return_codes[index] = command.wait()
                commands.remove((index, command))
                sleep_time = 0.05

        time.sleep(sleep_time)
        sleep_time = min(1, sleep_time * 2)

    if any(return_codes):
        out.write("Errors occured during processing!\n")
        out.flush()

    return return_codes
