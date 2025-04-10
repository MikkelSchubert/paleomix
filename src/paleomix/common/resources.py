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

import os
import shutil
import sys
from contextlib import AbstractContextManager
from importlib import resources
from pathlib import Path

from paleomix.common.fileutils import copy_file

if sys.version_info < (3, 11):
    from importlib.abc import Traversable
else:
    from importlib.resources.abc import Traversable


def _as_traversable(resource: str) -> Traversable:
    return resources.files("paleomix").joinpath("resources").joinpath(resource)


def read_text(resource: str) -> str:
    return _as_traversable(resource).read_text(encoding="utf-8")


def access(resource: str) -> AbstractContextManager[Path]:
    return resources.as_file(_as_traversable(resource))


def copy_resource(resource: str, destination: str) -> None:
    with access(resource) as path:
        if path.is_dir():
            shutil.copytree(path, destination)
        else:
            copy_file(path, destination)


def read_template(filename: str) -> str:
    """Returns the path to a report-file for a given tool."""
    return read_text(os.path.join("templates", filename))
