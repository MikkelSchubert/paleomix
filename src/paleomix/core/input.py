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
import sys
import termios
import tty
from typing import TYPE_CHECKING, Iterator

if TYPE_CHECKING:
    from typing_extensions import Self

    from paleomix.core.workers import HandleType

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

    def teardown(self) -> None:
        if self._tty_settings is not None:
            self._log.debug("Restoring terminal settings (enabling echo, etc.)")
            termios.tcsetattr(sys.stdin, termios.TCSADRAIN, self._tty_settings)
            self._tty_settings = None

    @property
    def handles(self) -> Iterator[HandleType]:
        if self._tty_settings is not None:
            yield sys.stdin.fileno()

    def process_key_presses(self) -> Iterator[CLIEvent]:
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
