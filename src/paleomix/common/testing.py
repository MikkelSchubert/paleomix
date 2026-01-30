# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import os
from pathlib import Path


class SetWorkingDirectory:
    """Sets the current working directory upon entry to that specified,
    in the constructor upon entry, and reverts to the previously used
    directory upon exiting a with statement."""

    def __init__(self, path: str | Path) -> None:
        self._old_cwd = None
        self._new_cwd = os.fspath(path)

    def __enter__(self) -> None:
        self._old_cwd = os.getcwd()
        os.chdir(self._new_cwd)

    def __exit__(self, typ: object, exc: object, tb: object) -> None:
        if self._old_cwd:
            os.chdir(self._old_cwd)
