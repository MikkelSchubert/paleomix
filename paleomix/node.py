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

import errno
import fnmatch
import itertools
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import paleomix
from paleomix.common import fileutils
from paleomix.common.command import AtomicCmd, CmdError, ParallelCmds, SequentialCmds
from paleomix.common.procs import quote_args
from paleomix.common.utilities import safe_coerce_to_frozenset
from paleomix.common.versions import Requirement

if TYPE_CHECKING:
    from paleomix.common.fileutils import PathTypes

_GLOBAL_ID = itertools.count()


class NodeError(RuntimeError):
    def __init__(self, *args: object, path: str | None = None) -> None:
        super().__init__(*args)
        self.path = path


class NodeMissingFilesError(NodeError):
    pass


class CmdNodeError(NodeError):
    pass


class NodeUnhandledError(NodeError):
    """This exception is thrown by Node.run() if a non-NodeError exception
    is raised in a subfunction (e.g. _setup, _run, or _teardown). The text
    for this exception will include both the original error message and a
    stacktrace for that error."""


class Node:
    __description: str | None
    input_files: frozenset[str]
    output_files: frozenset[str]
    intermediate_output_files: set[str]
    executables: frozenset[str]
    auxiliary_files: frozenset[str]
    requirements: frozenset[Requirement]

    threads: int
    dependencies: frozenset[Node]

    def __init__(
        self,
        description: str | None = None,
        threads: int = 1,
        input_files: Iterable[str] = (),
        output_files: Iterable[str] = (),
        executables: Iterable[str] = (),
        auxiliary_files: Iterable[str] = (),
        requirements: Iterable[Requirement] = (),
        dependencies: Iterable[Node] = (),
    ) -> None:
        if not (description is None or isinstance(description, str)):
            raise TypeError(description)

        self.__description = description
        self.input_files = self._validate_files(input_files)
        self.output_files = self._validate_files(output_files)
        self.intermediate_output_files = set()
        self.executables = self._validate_files(executables)
        self.auxiliary_files = self._validate_files(auxiliary_files)
        self.requirements = self._validate_requirements(requirements)

        self.threads = self._validate_nthreads(threads)
        self.dependencies = self._collect_nodes(dependencies)

        # If there are no input files, the node cannot be re-run based on
        # changes to the input, and nodes with output but no input are not
        # expected based on current usage.
        if not self.input_files and self.output_files:
            raise NodeError("Task not dependent upon input files: %s" % self)

        # Globally unique node ID
        self.id = next(_GLOBAL_ID)

    def run(self, temp_root: PathTypes) -> None:
        """Runs the node, by calling _setup, _run, and _teardown in that order.
        Prior to calling these functions, a temporary dir is created using the
        'temp_root' prefix from the config object. Both the config object and
        the temporary dir are passed to the above functions. The temporary
        dir is removed after _teardown is called, and all expected files
        should have been removed/renamed at that point.

        Any non-NodeError exception raised in this function is wrapped in a
        NodeUnhandledError, which includes a full backtrace. This is needed
        to allow showing these in the main process."""

        temp = None
        try:
            # Generate directory name and create dir at temp_root
            temp = self._create_temp_dir(temp_root)

            self._setup(temp)
            self._run(temp)
            self._teardown(temp)
            self._remove_temp_dir(temp)
        except NodeMissingFilesError:
            try:
                # The folder is most likely empty, but it is possible to re-use temp
                # directories for resumable tasks so we cannot delete it outrigth
                if temp is not None:
                    os.rmdir(temp)
            except OSError:
                pass

            raise
        except NodeError as error:
            self._write_error_log(temp, error)
            raise NodeError(
                "Error while running {}:\n  {}".format(
                    self, "\n  ".join(str(error).split("\n"))
                ),
                path=temp,
            ) from None
        except Exception as error:  # noqa: BLE001
            self._write_error_log(temp, error)
            raise NodeUnhandledError(
                f"Error while running {self}", path=temp
            ) from error

    def mark_intermediate_files(self, glob: str = "*") -> None:
        if glob == "*":
            self.intermediate_output_files = set(self.output_files)
            return

        self.intermediate_output_files.update(fnmatch.filter(self.output_files, glob))

    def _create_temp_dir(self, temp_root: PathTypes) -> str:
        """Called by 'run' in order to create a temporary folder."""
        return fileutils.create_temp_dir(temp_root)

    def _remove_temp_dir(self, temp: PathTypes) -> None:
        """Called by 'run' in order to remove an (now) empty temporary folder."""
        temp = fileutils.fspath(temp)
        log = logging.getLogger(__name__)
        for filename in self._collect_files(temp):
            log.warning(
                "Unexpected file in temporary directory: %r",
                os.path.join(temp, filename),
            )

        try:
            shutil.rmtree(temp)
        except OSError as error:
            if error.errno != errno.EBUSY:
                raise

            log.warning("Could not remove temporary directory: %r", error)

    def _setup(self, _temp: PathTypes) -> None:
        """Is called prior to '_run()' by 'run()'. Any code used to copy/link files,
        or other steps needed to ready the node for running may be carried out in this
        function. Checks that required input files exist, and raises an NodeError if
        this is not the case."""
        executables: list[str] = []
        for executable in self.executables:
            if executable == "%(PYTHON)s":
                executable = sys.executable

            executables.append(executable)

        missing_executables = fileutils.missing_executables(executables)
        if missing_executables:
            raise NodeError(f"Executable(s) not found: {missing_executables}")

        self._check_for_input_files(self.input_files | self.auxiliary_files)

    def _run(self, _temp: PathTypes) -> None:
        pass

    def _teardown(self, _temp: PathTypes) -> None:
        self._check_for_missing_files(self.output_files, "output")

    def __str__(self) -> str:
        """Returns the description passed to the constructor, or a default
        description if no description was passed to the constructor."""
        if self.__description:
            return self.__description
        return repr(self)

    def __getstate__(self) -> dict[str, Any]:
        """Called by pickle/cPickle to determine what to pickle; this is
        overridden to avoid pickling of requirements, dependencies, which would
        otherwise greatly inflate the amount of information that needs to be
        pickled."""
        obj_dict = self.__dict__.copy()
        obj_dict["requirements"] = ()
        obj_dict["dependencies"] = ()
        return obj_dict

    def _write_error_log(self, temp: str | None, error: Exception) -> None:
        if not (temp and os.path.isdir(temp)):
            return

        def _fmt(values: Iterable[str]) -> str:
            return "\n                   ".join(sorted(values))

        message = [
            f"PALEOMIX         = v{paleomix.__version__}",
            f"Command          = {quote_args(sys.argv)}",
            f"CWD              = {os.getcwd()!r}",
            "PATH             = {!r}".format(os.environ.get("PATH", "")),
            f"Task             = {self}",
            f"Threads          = {self.threads}",
            f"Input files      = {_fmt(self.input_files)}",
            f"Output files     = {_fmt(self.output_files)}",
            f"Auxiliary files  = {_fmt(self.auxiliary_files)}",
            f"Executables      = {_fmt(self.executables)}",
            "",
            f"Errors =\n{error}\n",
        ]
        message = "\n".join(message)

        try:
            with open(os.path.join(temp, "pipe.errors"), "w") as handle:
                handle.write(message)
        except OSError as oserror:
            sys.stderr.write(f"ERROR: Could not write failure log: {oserror}\n")

    def _collect_nodes(self, nodes: Iterable[Node]) -> frozenset[Node]:
        nodes = safe_coerce_to_frozenset(nodes)
        for node in nodes:
            if not isinstance(node, Node):
                raise TypeError(node)

        return nodes

    def _check_for_input_files(self, filenames: Iterable[str]) -> None:
        missing_files = fileutils.missing_files(filenames)
        if missing_files:
            raise NodeMissingFilesError(
                "Missing input files for command:\n"
                "\t- Command: {}\n"
                "\t- Files: {}".format(self, "\n\t         ".join(missing_files))
            )

    def _check_for_missing_files(
        self,
        filenames: Iterable[str],
        description: str,
    ) -> None:
        missing_files = fileutils.missing_files(filenames)
        if missing_files:
            message = (
                "Missing {} files for command:\n\t- Command: {}\n\t- Files: {}".format(
                    description, self, "\n\t         ".join(missing_files)
                )
            )
            raise NodeError(message)

    @classmethod
    def _validate_requirements(
        cls, requirements: Iterable[Requirement]
    ) -> frozenset[Requirement]:
        requirements = safe_coerce_to_frozenset(requirements)
        for requirement in requirements:
            if not isinstance(requirement, Requirement):
                raise TypeError(requirement)
        return requirements

    @classmethod
    def _validate_files(cls, files: Iterable[str]) -> frozenset[str]:
        return frozenset(fileutils.validate_filenames(files))

    @classmethod
    def _validate_nthreads(cls, threads: object) -> int:
        if not isinstance(threads, int):
            raise TypeError(f"'threads' must be a positive integer, not {threads!r}")
        elif threads < 1:
            raise ValueError(
                "'threads' must be a positive integer, not %i" % (threads,)
            )
        return threads

    @staticmethod
    def _collect_files(root: PathTypes) -> Iterable[str]:
        root = fileutils.fspath(root)

        def _walk_dir(path: str) -> Iterable[str]:
            for entry in os.scandir(path):
                if entry.is_file():
                    yield str(Path(entry.path).relative_to(root))
                elif entry.is_dir():
                    yield from _walk_dir(entry.path)

        yield from _walk_dir(root)


class CommandNode(Node):
    _command: AtomicCmd | ParallelCmds | SequentialCmds

    def __init__(
        self,
        command: AtomicCmd | ParallelCmds | SequentialCmds,
        description: str | None = None,
        threads: int = 1,
        dependencies: Iterable[Node] = (),
    ) -> None:
        Node.__init__(
            self,
            description=description,
            input_files=command.input_files,
            output_files=command.output_files,
            auxiliary_files=command.auxiliary_files,
            executables=command.executables,
            requirements=command.requirements,
            threads=threads,
            dependencies=dependencies,
        )

        self._command = command

    def _run(self, temp: PathTypes) -> None:
        """Runs the command object provided in the constructor, and waits for it to
        terminate. If any errors during the running of the command, this function
        raises a NodeError detailing the returned error-codes."""
        try:
            self._command.run(temp)
        except CmdError as error:
            raise CmdNodeError(f"{self._command!s}\n\n{error}") from error

        return_codes = self._command.join()
        if any(return_codes):
            raise CmdNodeError(str(self._command))

    def _teardown(self, temp: PathTypes) -> None:
        required_files = self._command.expected_temp_files
        current_files = set(self._collect_files(temp))

        missing_files = required_files - current_files
        if missing_files:
            raise CmdNodeError(
                (
                    "Error running task, required files were not created:\n"
                    "Temporary directory: {!r}\n"
                    "\tRequired files missing from temporary directory:\n\t    - {}"
                ).format(temp, "\n\t    - ".join(sorted(map(repr, missing_files))))
            )

        self._command.commit()

        super()._teardown(temp)
