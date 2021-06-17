#!/usr/bin/python
#
# Copyright (c) 2021 Mikkel Schubert <MikkelSch@gmail.com>
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
import logging
import os
import select
import sys
import termios
import tty
from typing import Any, Iterator, List

_COMMANDS = {
    "h": "Prints this message.",
    "l": "Lists the currently runnning tasks.",
    "+": "Increases the maximum number of threads by one.",
    "-": "Decreases the maximum number of threads by one; does not kill running tasks.",
}


class CLIEvent:
    pass


class ListTasksEvent(CLIEvent):
    pass


class ThreadsEvent(CLIEvent):
    def __init__(self, change: int):
        self.change = change


class CommandLine(object):
    def __init__(self):
        self._tty_settings = None
        self._log = logging.getLogger(__name__)

    def __enter__(self):
        self.setup()

        return self

    def __exit__(self, _type: Any, _value: Any, _traceback: Any):
        self.teardown()

    @property
    def active(self):
        return self._tty_settings is not None

    def setup(self):
        if self._tty_settings is None:
            # False if the pipeline is being piped somewhere
            if sys.stdin.isatty() and sys.stdout.isatty():
                # False if the process is running in the background
                if os.getpgrp() == os.tcgetpgrp(sys.stdout.fileno()):
                    try:
                        # Store old settings
                        self._tty_settings = termios.tcgetattr(sys.stdin)
                        # Disable echo
                        tty.setcbreak(sys.stdin.fileno())
                    except termios.error:
                        pass  # Silently ignore failures

    def teardown(self):
        if self._tty_settings is not None:
            # Restore settings (re-enable echo)
            termios.tcsetattr(sys.stdin, termios.TCSADRAIN, self._tty_settings)
            self._tty_settings = None

    @property
    def handles(self) -> List[Any]:
        if self._tty_settings is None:
            return []

        return [sys.stdin]

    def process_key_presses(self) -> Iterator[CLIEvent]:
        while self._tty_settings and self._poll_stdin():
            character = sys.stdin.read(1)

            if character == "+":
                yield ThreadsEvent(1)
            elif character == "-":
                yield ThreadsEvent(-1)
            elif character in "lL":
                yield ListTasksEvent()
            elif character in "hH":
                self._log_help()

    @classmethod
    def _poll_stdin(cls) -> bool:
        return select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], [])

    def _log_help(self):
        self._log.info("Commands:")
        self._log.info("  Key   Function")
        for key, help in _COMMANDS.items():
            self._log.info("  %s    %s", key, help)
