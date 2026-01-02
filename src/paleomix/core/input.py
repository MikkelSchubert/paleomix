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

import logging
import os
import select
import signal
import socket
import sys
import termios
import tty
from collections.abc import Iterator
from multiprocessing.connection import Connection
from typing import Union

from typing_extensions import Self, TypeAlias

HandleType: TypeAlias = Union[Connection, socket.socket, int]


_COMMANDS = {
    "h": "Prints this message.",
    "l": "Lists the currently running tasks.",
    "+": "Increases the maximum number of threads by one.",
    "-": "Decreases the maximum number of threads by one; does not kill running tasks.",
    "Ctrl+C": "Quit after finishing running tasks; press twice to quit immediately.",
}


class CLIEvent:
    pass


class ListTasksEvent(CLIEvent):
    pass


class ThreadsEvent(CLIEvent):
    def __init__(self, change: int) -> None:
        self.change = change


class CommandLine:
    def __init__(self) -> None:
        self._tty_settings = None
        self._log = logging.getLogger(__name__)

        # Self-socket used to interrupt polling on interrupts
        self._csock = self._csock = None

    def __enter__(self) -> Self:
        self.setup()

        return self

    def __exit__(self, typ: object, exc: object, tb: object) -> None:
        self.teardown()

    @property
    def active(self) -> bool:
        return self._tty_settings is not None

    def setup(self) -> None:
        if (
            self._tty_settings is None
            # False if the pipeline is being piped somewhere
            and sys.stdin.isatty()
            and sys.stdout.isatty()
            # False if the process is running in the background
            and os.getpgrp() == os.tcgetpgrp(sys.stdout.fileno())
        ):
            try:
                # Store old settings
                self._tty_settings = termios.tcgetattr(sys.stdin)
                # Disable echo
                tty.setcbreak(sys.stdin.fileno())
            except termios.error:
                self._tty_settings = None
                self._log.debug("Command-line interface could not be enabled")
            else:
                self._log.debug("Enabled command-line interface")
        else:
            self._log.debug("Command-line interface disabled; not in interactive mode")

        if self._csock is None and self._csock is None:
            self._ssock, self._csock = socket.socketpair()
            self._ssock.setblocking(False)
            self._csock.setblocking(False)

            signal.set_wakeup_fd(self._csock.fileno(), warn_on_full_buffer=False)

    def teardown(self) -> None:
        if self._tty_settings is not None:
            self._log.debug("Restoring terminal settings (enabling echo, etc.)")
            termios.tcsetattr(sys.stdin, termios.TCSADRAIN, self._tty_settings)
            self._tty_settings = None

        if self._ssock is not None:
            signal.set_wakeup_fd(-1)
            self._ssock.close()
            self._ssock = None

        if self._csock is not None:
            self._csock.close()
            self._csock = None

    @property
    def handles(self) -> Iterator[HandleType]:
        if self._tty_settings is not None:
            yield sys.stdin.fileno()

        if self._ssock is not None:
            yield self._ssock.fileno()

    def process_key_presses(self) -> Iterator[CLIEvent]:
        while self._ssock is not None:
            try:
                data = self._ssock.recv(8)
                if not data:
                    break
            except InterruptedError:
                continue
            except BlockingIOError:
                break

        if self._tty_settings is not None:
            characters: list[str] = []
            while self._poll_stdin():
                characters.append(sys.stdin.read(1))

            character = "".join(characters)
            if character == "+":
                yield ThreadsEvent(1)
            elif character == "-":
                yield ThreadsEvent(-1)
            elif character in ["l", "L"]:
                yield ListTasksEvent()
            elif character in ["h", "H"]:
                self._log_help()

    @classmethod
    def _poll_stdin(cls) -> bool:
        return select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], [])

    def _log_help(self) -> None:
        self._log.info("Commands:")
        self._log.info("  Key      Function")
        for key, text in _COMMANDS.items():
            self._log.info("  %s  %s", key.ljust(7), text)
