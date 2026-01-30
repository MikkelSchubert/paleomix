# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
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
