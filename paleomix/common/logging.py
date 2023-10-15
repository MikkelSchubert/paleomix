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

import copy
import errno
import itertools
import logging
import os
import sys
import time
from logging import LogRecord
from typing import TYPE_CHECKING, Iterator

import coloredlogs
from humanfriendly.terminal import ansi_wrap, terminal_supports_colors

if TYPE_CHECKING:
    from io import TextIOWrapper

    from paleomix.common.argparse import ArgumentParser


_CONSOLE_MESSAGE_FORMAT: str = "%(asctime)s %(levelname)s %(status)s%(message)s"
_FILE_MESSAGE_FORMAT: str = "%(asctime)s %(name)s %(levelname)s %(status)s%(message)s"


class LogRecordWithStatus(logging.LogRecord):
    status: str | Status | None


class BasicFormatter(coloredlogs.ColoredFormatter):
    def format(self, record: LogRecord) -> str:  # noqa: A003
        record = copy.copy(record)
        record.status = self._process_status(record)

        message = record.msg
        if not isinstance(message, str):
            message = str(message)

        if record.args:
            message = message % record.args

        if "\n" in message:
            lines: list[str] = []
            record.args = ()
            for line in message.split("\n"):
                record.msg = line
                lines.append(super().format(record))

            return "\n".join(lines)

        return super().format(record)

    def _process_status(self, record: LogRecord) -> str:
        if isinstance(record, LogRecordWithStatus):
            return f"[{record.status}] "

        return ""


class PaleomixFormatter(BasicFormatter):
    def _process_status(self, record: LogRecord) -> str:
        status = getattr(record, "status", None)
        if isinstance(status, Status) and status.color:
            return f"[{ansi_wrap(str(status), color=status.color)}] "

        return ""


def initialize_console_logging(log_level: str = "info") -> None:
    log_level = coloredlogs.level_to_number(log_level)  # type: ignore

    logger = logging.getLogger()
    logger.setLevel(log_level)

    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler) and handler.stream is sys.stderr:
            break
    else:
        handler = logging.StreamHandler()
        logger.addHandler(handler)

    fmt_class = PaleomixFormatter if terminal_supports_colors() else BasicFormatter
    formatter = fmt_class(fmt=_CONSOLE_MESSAGE_FORMAT)

    handler.setFormatter(formatter)


def initialize(
    log_level: str = "info",
    log_file: str | None = None,
    auto_log_file: str | None = "paleomix",
) -> None:
    initialize_console_logging(log_level)

    if log_file:
        logger = logging.getLogger(__name__)
        logger.info("Writing %s log to %r", log_level, log_file)

        handler = logging.FileHandler(log_file)
        handler.setFormatter(BasicFormatter(_FILE_MESSAGE_FORMAT))
        handler.setLevel(coloredlogs.level_to_number(log_level))

        root = logging.getLogger()
        root.addHandler(handler)
    elif auto_log_file:
        template = "{}.{}_%02i.log".format(
            auto_log_file, time.strftime("%Y%m%d_%H%M%S")
        )
        handler = LazyLogfile(template, log_level=logging.ERROR)
        handler.setFormatter(BasicFormatter(_FILE_MESSAGE_FORMAT))
        handler.setLevel(logging.ERROR)

        root = logging.getLogger()
        root.addHandler(handler)


def add_argument_group(parser: ArgumentParser, *, log_file: bool = True) -> None:
    """Adds an option-group to an OptionParser object, with options
    pertaining to logging. Note that 'initialize' expects the config
    object to have these options."""
    group = parser.add_argument_group("Logging")
    if log_file:
        group.add_argument(
            "--log-file",
            default=None,
            help="Write log-messages to this file. If a log-file is not specified, "
            "error messages, and only error messages, will be written to an "
            "automatically generated log file",
        )

    group.add_argument(
        "--log-level",
        default="info",
        choices=("debug", "info", "warning", "error"),
        type=str.lower,
        help="Log messages at the specified level. This option applies to the "
        "`--log-file` option and to log messages printed to the terminal.",
    )


def get_logfiles() -> Iterator[str]:
    root = logging.getLogger()
    for handler in root.handlers:
        if isinstance(handler, logging.FileHandler) and handler.stream:
            yield handler.baseFilename


class LazyLogfile(logging.FileHandler):
    def __init__(self, template: str, log_level: int) -> None:
        logging.FileHandler.__init__(self, template, delay=True)
        # Use absolute path for template
        self._template = self.baseFilename
        self._log_level = log_level

    def emit(self, record: logging.LogRecord) -> None:
        # Don't try to log self-emitted log events, to avoid recursive loops
        if record.name != __name__:
            super().emit(record)

    def _open(self) -> TextIOWrapper:
        """Try to open a new logfile, taking steps to ensure that
        existing logfiles using the same template are not clobbered."""
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL

        for start in itertools.count(start=1):
            filename = self._template % (start,)
            os.makedirs(os.path.dirname(filename), exist_ok=True)

            try:
                stream = os.fdopen(os.open(filename, flags), "w")
                logger = logging.getLogger(__name__)
                level_name = logging.getLevelName(self._log_level).lower()
                logger.info("Saving %s logs to %r", level_name, filename)

                self.baseFilename = filename
            except OSError as error:
                if error.errno != errno.EEXIST:
                    raise
            else:
                return stream

        raise AssertionError("loop should have ended")


class Status:
    def __init__(self, color: str | None = None) -> None:
        self.color = color

    def __str__(self) -> str:
        raise NotImplementedError
