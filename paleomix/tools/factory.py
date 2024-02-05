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
"""Factory for AtomicCmds for the various PALEOMIX commands.

Ensures that the version called corresponds to the running version, in case
multiple versions are present in the users' PATH, or that the current version
is not available from the users' PATH.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, TypeVar

from paleomix.common.command import (
    ArgsType,
    AtomicCmd,
    AtomicFileTypes,
    InputFile,
    OutputFile,
)
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.versions import Requirement

if TYPE_CHECKING:
    from pathlib import Path

PYTHON_VERSION = Requirement(
    ["%(PYTHON)s", "--version"],
    regexp=r"Python (\d+\.\d+\.\d+)",
    specifiers=">=3.7",
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
