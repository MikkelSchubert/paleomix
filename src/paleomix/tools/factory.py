# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
"""Factory for AtomicCmds for the various PALEOMIX commands.

Ensures that the version called corresponds to the running version, in case
multiple versions are present in the users' PATH, or that the current version
is not available from the users' PATH.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import TypeVar

from paleomix.common.command import (
    ArgsType,
    AtomicCmd,
    AtomicFileTypes,
    InputFile,
    OutputFile,
)
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.versions import Requirement

PYTHON_VERSION = Requirement(
    ["%(PYTHON)s", "--version"],
    regexp=r"Python (\d+\.\d+\.\d+)",
    specifiers=">=3.8",
    name="Python",
)

RSCRIPT_VERSION = Requirement(
    call=("Rscript", "--version"),
    regexp=r"version (\d+\.\d+\.\d+)",
    specifiers=">=3.3.3",
)

# Avoid returning a union type in the the case where all values are str
T = TypeVar("T", str, ArgsType)


def command(args: Iterable[T]) -> tuple[T | str, ...]:
    return ("%(PYTHON)s", "-m", "paleomix", *safe_coerce_to_tuple(args))


def new(
    args: Iterable[ArgsType],
    *,
    stdin: int | str | Path | InputFile | AtomicCmd | None = None,
    stdout: int | str | Path | OutputFile | None = None,
    stderr: int | str | Path | OutputFile | None = None,
    set_cwd: bool = False,
    extra_files: Iterable[AtomicFileTypes] = (),
    requirements: Iterable[Requirement] = (),
) -> AtomicCmd:
    """Returns AtomicCmd setup to call the tools accessible through the
    'paleomix' command-line tool. This builder adds executable / version checks
    for the specified command, but does not add any arguments. Thus, calling
    new with the argument "cat" produces the equivalent of ["paleomix", "cat"].
    """
    return AtomicCmd(
        command(args),
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        set_cwd=set_cwd,
        extra_files=extra_files,
        requirements=[*requirements, PYTHON_VERSION],
    )


def rscript_command(args: Iterable[T] = ()) -> tuple[T | str, ...]:
    return command((":rscript", *safe_coerce_to_tuple(args)))


def rscript(
    args: Iterable[ArgsType] = (),
    *,
    stdin: int | str | Path | InputFile | AtomicCmd | None = None,
    stdout: int | str | Path | OutputFile | None = None,
    stderr: int | str | Path | OutputFile | None = None,
    set_cwd: bool = False,
    extra_files: Iterable[AtomicFileTypes] = (),
    requirements: Iterable[Requirement] = (),
) -> AtomicCmd:
    return new(
        (":rscript", *safe_coerce_to_tuple(args)),
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        set_cwd=set_cwd,
        extra_files=extra_files,
        requirements=[*requirements, RSCRIPT_VERSION],
    )
