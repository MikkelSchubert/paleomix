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
import multiprocessing
import os
import select
import sys
import termios
import tty


_COMMANDS = {
    "h": "Prints this message.",
    "l": "Lists the currently runnning tasks.",
    "+": "Increases the maximum number of threads by one.",
    "-": "Decreases the maximum number of threads by one; does not kill running tasks.",
}


class CommandLine(object):
    def __init__(self):
        self._tty_settings = None
        self._log = logging.getLogger(__name__)

    def __enter__(self):
        self.setup()

        return self

    def __exit__(self, _type, _value, _traceback):
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
                    except tty.error:
                        pass  # Silently ignore failures

    def teardown(self):
        if self._tty_settings is not None:
            # Restore settings (re-enable echo)
            termios.tcsetattr(sys.stdin, termios.TCSADRAIN, self._tty_settings)
            self._tty_settings = None

    @property
    def handles(self):
        if self._tty_settings is None:
            return []

        return [sys.stdin]

    def process_key_presses(self, threads=1, workers=()):
        if not self._tty_settings:
            return threads

        old_threads = threads
        while self.poll_stdin():
            character = sys.stdin.read(1)
            if character == "+":
                threads = min(multiprocessing.cpu_count(), threads + 1)
            elif character == "-":
                threads = max(0, threads - 1)
            elif character in "lL":
                self._log_tasks(workers)
            elif character in "hH":
                self._log.info("Commands:")
                self._log.info("  Key   Function")
                for key, help in _COMMANDS.items():
                    self._log.info("  %s    %s", key, help)

        if threads != old_threads:
            self._log.info("Max threads changed from %i to %i", old_threads, threads)

        return threads

    @classmethod
    def poll_stdin(cls):
        return select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], [])

    def _log_tasks(self, workers):
        total_threads = 0
        total_workers = 0
        for worker in workers:
            tasks = sorted(worker.tasks, key=lambda it: it.id)
            threads = sum(task.threads for task in tasks)

            if tasks:
                self._log.info(
                    "Running %i tasks on %s (using %i/%i threads):",
                    len(tasks),
                    worker.name,
                    threads,
                    worker.threads,
                )

                for idx, task in enumerate(tasks, start=1):
                    self._log.info("  % 2i. %s", idx, task)
            else:
                self._log.info(
                    "No tasks running on %s (using 0/%i threads)",
                    worker.name,
                    worker.threads,
                )

            total_threads += threads
            total_workers += 1

        if total_workers > 1:
            self._log.info(
                "A total of %i threads are used across %i workers",
                total_threads,
                total_workers,
            )
