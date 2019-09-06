#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MikkelSch@gmail.com>
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
import errno
import itertools
import logging
import optparse
import os
import time

import coloredlogs


_LOG_LEVELS = {
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "debug": logging.DEBUG,
}

_LOG_FORMAT = "%(asctime)s %(name)s %(levelname)s %(message)s"


def initialize(log_level="info", log_file=None, name="paleomix"):
    log_level = _LOG_LEVELS[log_level.lower()]

    coloredlogs.install(fmt=_LOG_FORMAT)
    root = logging.getLogger()

    if log_file:
        handler = logging.FileHandler(log_file)
        handler.setFormatter(logging.Formatter(_LOG_FORMAT))
        root.addHandler(handler)
    elif name:
        template = '%s.%s_%%02i.log' % (name, time.strftime("%Y%m%d_%H%M%S"))
        handler = LazyLogfile(template)
        handler.setFormatter(logging.Formatter(_LOG_FORMAT))
        handler.setLevel(logging.WARNING)
        root.addHandler(handler)


def add_optiongroup(parser, default="info"):
    """Adds an option-group to an OptionParser object, with options
    pertaining to logging. Note that 'initialize' expects the config
    object to have these options."""
    group = optparse.OptionGroup(parser, "Logging")
    group.add_option(
        "--log-file",
        default=None,
        help="Create the specified log-file and write any messages to this file.",
    )
    group.add_option(
        "--log-level",
        default=default,
        type="choice",
        choices=("info", "warning", "error", "debug"),
        help="Log messages to log-file at and above the specified level; one of "
        "'info', 'warning', 'error', or 'debug' [%default]",
    )
    parser.add_option_group(group)


def get_logfiles():
    root = logging.getLogger()
    for handler in root.handlers:
        if isinstance(handler, logging.FileHandler) and handler.stream:
            yield handler.baseFilename


class LazyLogfile(logging.FileHandler):
    def __init__(self, template):
        logging.FileHandler.__init__(self, template, delay=True)
        # Use absolute path for template
        self._template = self.baseFilename

    def _open(self):
        """Try to open a new logfile, taking steps to ensure that
        existing logfiles using the same template are not clobbered."""
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL

        for start in itertools.count(start=1):
            filename = self._template % (start,)

            try:
                stream = os.fdopen(os.open(filename, flags), "w")
                logger = logging.getLogger(__name__)
                logger.info('Logging warnings and errors to %r', filename)

                self.baseFilename = filename
                return stream
            except OSError as error:
                if error.errno != errno.EEXIST:
                    raise
