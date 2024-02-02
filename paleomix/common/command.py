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

import collections
import os
import shlex
import signal
import subprocess
import sys
from pathlib import Path
from typing import IO, TYPE_CHECKING, Callable, Dict, Iterable, List, Tuple, Union

from paleomix.common import fileutils
from paleomix.common.procs import RegisteredPopen
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.versions import Requirement

if TYPE_CHECKING:
    from typing_extensions import Self

    from paleomix.common.fileutils import PathTypes


class CmdError(RuntimeError):
    """Exception raised for AtomicCmd specific errors."""

    def __init__(self, msg: object) -> None:
        RuntimeError.__init__(self, msg)


class _AtomicFile:
    def __init__(self, path: PathTypes) -> None:
        self.path = fileutils.fspath(path)

        if isinstance(self.path, bytes):
            raise TypeError(f"invalid path {path!r}")

    def basename(self) -> str:
        return os.path.basename(self.path)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.path!r})"


class AuxiliaryFile(_AtomicFile):
    pass


class Executable(_AtomicFile):
    pass


class _IOFile(_AtomicFile):
    def __init__(self, path: PathTypes, *, temporary: bool = False) -> None:
        super().__init__(path)
        self.temporary = bool(temporary)

        if temporary and os.path.dirname(self.path):
            raise ValueError(f"directory component in temporary path {self.path!r}")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.path!r}, {self.temporary})"


class InputFile(_IOFile):
    pass


class TempInputFile(InputFile):
    def __init__(self, path: PathTypes) -> None:
        super().__init__(os.path.basename(path), temporary=True)


class OutputFile(_IOFile):
    pass


class TempOutputFile(OutputFile):
    def __init__(self, path: PathTypes) -> None:
        super().__init__(os.path.basename(path), temporary=True)


IOFileTypes = Union[InputFile, OutputFile, TempInputFile, TempOutputFile]
AtomicFileTypes = Union[AuxiliaryFile, Executable, IOFileTypes]

# Possible types of .stdin/.stdout/.stderr, int being either DEVNULL or PIPE
WrappedPipeType = Union[int, IOFileTypes, "AtomicCmd"]
# Types that can be passed as values for STDIN, STDOUT, and STDERR
PipeType = Union[None, str, Path, WrappedPipeType]

# Pos
OptionValueType = Union[str, float, IOFileTypes, None]
OptionsType = Dict[
    str,
    Union[
        OptionValueType,
        List[OptionValueType],
        Tuple[OptionValueType, ...],
    ],
]

JoinType = List[Union[str, None, int]]


class AtomicCmd:
    """Executes a command, only moving resulting files to the destination directory if
    the command was successful. This prevent the accidental use of partial files
    in downstream analysis and eases restarting pipelines following errors (no cleanup).

    Individual files/paths in the command are specified using the InputFile, OutputFile,
    AuxiliaryFile, and Executable classes, allowing the command to be transparently
    modified to execute in a temporary directory.

    When an AtomicCmd is run(), a signal handler is installed for SIGTERM, which
    ensures that any running processes are terminated. In the absence of this,
    AtomicCmds run in terminated subprocesses can result in still running children
    after the termination of the parents."""

    PIPE = subprocess.PIPE
    DEVNULL = subprocess.DEVNULL

    _command: list[str | AtomicFileTypes]
    _proc: subprocess.Popen[bytes] | None
    _temp: str | None
    _running: bool
    _set_cwd: bool
    _terminated: bool

    _executables: set[str]
    _requirements: set[Requirement]
    _input_files: set[IOFileTypes]
    _output_files: set[IOFileTypes]
    _output_basenames: set[str]
    _auxiliary_files: set[str]

    _stdin: WrappedPipeType
    _stdout: WrappedPipeType
    _stderr: WrappedPipeType

    def __init__(
        self,
        command: Iterable[str | int | Path | AtomicFileTypes],
        *,
        stdin: int | str | Path | InputFile | AtomicCmd | None = None,
        stdout: int | str | Path | OutputFile | None = None,
        stderr: int | str | Path | OutputFile | None = None,
        set_cwd: bool = False,
        extra_files: Iterable[AtomicFileTypes] = (),
        requirements: Iterable[Requirement] = (),
    ) -> None:
        """Takes a command and a set of files.

        The command is expected to be an iterable starting with the name of an
        executable, with each item representing one string on the command line. Thus,
        the command "find /etc -name 'profile*'" might be represented as the list
        ["find", "/etc", "-name", "profile*"].

        Commands typically consist of an executable, one or more input files, one or
        more output files, and zero or more pipes. Files are specified using the
        InputFile, OutputFile, AuxiliaryFile, and Executable classes, which allows easy
        tracking of requirements and other features. Note that only files, and not
        directories, are supported as input/output!

        Keyword arguments:
            stdin -- A string path, InputFile, DEVNULL, or an AtomicCmd instance.
                     If stdin is an AtomicCmd instance, stdin is read from its stdout.
            stdout -- A string path, OutputFile, PIPE, or DEVNULL.
            stderr -- A string path, OutputFile, or DEVNULL.
            set_cwd -- If true, cwd is set to the temporary folder before executing the
                       command. InputFiles, OutputFiles, and AuxiliaryFiles are
                       automatically converted to paths relative to the temporary dir.
            extra_files -- An optional sequence of InputFiles, OutputFiles,
                           AuxiliaryFiles and Executables used by the command, that are
                           not explicitly part of the command line arguments. This
                           could include data files and executables indirectly executed
                           by a script.
            requirements -- An optional sequence of Requirements that can be used to
                            check that requirements are met for invoking this command.

        EXAMPLE 1: Creating a gzipped tar-archive from two files
        The command "tar cjf /path/to/output.tbz2 /path/to/input1 /path/to/input2"
        could be represented using the following AtomicCmd:
        cmd = AtomicCmd(["tar", "cjf", OutputFile("/path/to/output.tbz2"),
                         InputFile("/path/to/input1"), InputFile("/path/to/input2")])

        InputFiles and OutputFiles marked as 'temporary' are used for files that are
        read from or written to the temporary directory and which are deleted upon
        completion of the command (when calling 'commit'). Such files are not allowed
        to have directory component.

        EXAMPLE 2: zcat'ing an archive
        The command "zcat /path/to/input > output" could be represented as follows:
        cmd = AtomicCmd(["zcat", InputFile("/path/to/input")], stdout="output")
        """
        self._command = []
        self._proc = None
        self._temp = None
        self._running = False
        self._set_cwd = set_cwd
        self._terminated = False

        self._executables = set()
        self._requirements = set(requirements)
        self._input_files = set()
        self._output_files = set()
        self._output_basenames = set()
        self._auxiliary_files = set()

        self.append(*safe_coerce_to_tuple(command))
        if not self._command or not self._command[0]:
            raise ValueError("Empty command in AtomicCmd constructor")
        elif isinstance(self._command[0], str):
            self._executables.add(self._command[0])
            # Ensure that executables with path components are handled properly
            self._command[0] = Executable(self._command[0])
        elif not isinstance(self._command[0], Executable):
            raise TypeError(f"exe must be str or Executable, not {self._command[0]}")

        self._stdin = self._wrap_pipe(InputFile, stdin)
        self._stdout = self._wrap_pipe(OutputFile, stdout, "stdout")
        self._stderr = self._wrap_pipe(OutputFile, stderr, "stderr")

        self.add_extra_files(extra_files)
        self.add_extra_files(
            pipe
            for pipe in (self._stdin, self._stdout, self._stderr)
            if isinstance(pipe, _AtomicFile)
        )

        for value in self._requirements:
            if not isinstance(value, Requirement):
                raise TypeError(value)

    def append(self, *args: AtomicFileTypes | str | float) -> None:
        if self._proc is not None:
            raise CmdError("cannot modify already started command")

        for value in args:
            if value == "%(PYTHON)s":
                value = Executable("%(PYTHON)s")

            if isinstance(value, _AtomicFile):
                self._record_atomic_file(value)
            else:
                value = str(value)

            self._command.append(value)

    def add_requirement(self, requirements: Requirement) -> None:
        self._requirements.add(requirements)

    def add_extra_files(self, files: Iterable[AtomicFileTypes]) -> None:
        if self._proc is not None:
            raise CmdError("cannot modify already started command")

        for value in files:
            if not isinstance(value, _AtomicFile):
                raise TypeError(value)

            self._record_atomic_file(value)

    def append_options(
        self,
        options: OptionsType,
        pred: Callable[[str], bool] = lambda s: s.startswith("-"),
    ) -> None:
        if not isinstance(options, dict):
            raise TypeError(f"options must be dict, not {options!r}")

        for key, values in options.items():
            if not isinstance(key, str):
                raise TypeError(f"keys must be strings, not {key!r}")
            elif not pred(key):
                continue

            if isinstance(values, (list, tuple)):
                for value in values:
                    if not isinstance(value, (int, str, float, _AtomicFile)):
                        raise TypeError(value)

                    self.append(key)
                    self.append(value)
            elif values is None:
                self.append(key)
            elif isinstance(values, (int, str, float, _AtomicFile)):
                self.append(key)
                self.append(values)
            else:
                raise TypeError(values)

    def merge_options(
        self,
        *,
        user_options: OptionsType | None = None,
        fixed_options: OptionsType | None = None,
        blacklisted_options: Iterable[str] = (),
        pred: Callable[[str], bool] = lambda s: s.startswith("-"),
    ) -> None:
        if not isinstance(fixed_options, dict):
            raise TypeError(f"options must be dict, not {fixed_options!r}")
        elif not isinstance(user_options, dict):
            raise TypeError(f"user_options must be dict, not {user_options!r}")

        errors: list[str] = []
        for key in user_options.keys() & fixed_options.keys():
            errors.append(f"{key} cannot be overridden")

        for key in user_options.keys() & blacklisted_options:
            errors.append(f"{key} is not supported")

        if errors:
            raise CmdError(
                "invalid command-line options for {!r}: {}".format(
                    " ".join(self.to_call("%(TEMP_DIR)s")), "\n".join(errors)
                )
            )

        self.append_options(fixed_options, pred=pred)
        self.append_options(user_options, pred=pred)

    @property
    def output_files(self) -> set[str]:
        return {it.path for it in self._output_files if not it.temporary}

    @property
    def expected_temp_files(self) -> set[str]:
        return {it.basename() for it in self._output_files if not it.temporary}

    @property
    def optional_temp_files(self) -> set[str]:
        return {it.basename() for it in self._output_files if it.temporary}

    @property
    def input_files(self) -> set[str]:
        return {it.path for it in self._input_files if not it.temporary}

    @property
    def executables(self) -> set[str]:
        return set(self._executables)

    @property
    def auxiliary_files(self) -> set[str]:
        return set(self._auxiliary_files)

    @property
    def requirements(self) -> set[Requirement]:
        return set(self._requirements)

    def run(self, temp: PathTypes) -> None:
        """Runs the given command, saving files in the specified temp folder.
        To move files to their final destination, call commit(). Note that in
        contexts where the *Cmds classes are used, this function may block.

        """
        if self._running:
            raise CmdError("Calling 'run' on already running command.")

        temp = fileutils.fspath(temp)
        self._temp = temp
        self._running = True

        stdin = stdout = stderr = None
        try:
            stdin = self._open_pipe(temp, self._stdin, "rb")
            stdout = self._open_pipe(temp, self._stdout, "wb")
            stderr = self._open_pipe(temp, self._stderr, "wb")

            cwd = temp if self._set_cwd else None
            call = self.to_call(temp)

            self._proc = RegisteredPopen(
                call,
                stdin=stdin,
                stdout=stdout,
                stderr=stderr,
                cwd=cwd,
                start_new_session=True,
            )
        except (OSError, ValueError, CmdError) as error:
            message = "Error running commands:\n  Call = {!r}\n  Error = {!r}"
            raise CmdError(message.format(self._command, error)) from error
        finally:
            # Close pipes to allow the command to receive SIGPIPE
            for handle in (stdin, stdout, stderr):
                if not (handle is None or isinstance(handle, int)):
                    handle.close()

    def communicate(self) -> tuple[bytes, bytes]:
        if not self._running:
            raise CmdError("Calling 'communicate' on non-running command.")
        elif self._proc is None:
            raise AssertionError("self._proc is none in AtomicCmd.communicate")

        return self._proc.communicate()

    def ready(self) -> bool:
        """Returns true if the command has been run to completion,
        regardless of wether or not an error occured."""
        return bool(self._proc and self._proc.poll() is not None)

    def join(self, timeout: float | None = None) -> JoinType:
        """Similar to Popen.wait(), but returns the value wrapped in a list.
        Must be called before calling commit."""
        if not self._proc:
            return [None]

        return_code = self._proc.wait(timeout)
        self._running = False

        if return_code < 0:
            return_code = signal.Signals(-return_code).name
        return [return_code]

    def wait(self) -> str | int | None:
        """Equivalent to Subprocess.wait. This function should only
        be used in contexts where a AtomicCmd needs to be combined
        with Subprocesses, as it does not exist for AtomicSets."""
        return self.join()[0]

    def terminate(self, signal: int = signal.SIGTERM) -> None:
        """Sends SIGTERM to process if it is still running.
        Has no effect if the command has already finished."""
        if self._proc and self._proc.poll() is None:
            try:
                os.killpg(self._proc.pid, signal)
                self._terminated = True
            except OSError:
                pass  # Already dead / finished process

    def commit(self) -> None:
        if not self.ready():
            raise CmdError("Attempting to commit before command has completed")
        elif self._temp is None:
            raise AssertionError("self._temp should not be None")

        if self._running:
            raise CmdError("Called 'commit' before calling 'join'")

        expected_files = {os.path.basename(fpath) for fpath in self.output_files}
        missing_files = expected_files - set(os.listdir(self._temp))
        if missing_files:
            raise CmdError(
                "Expected files not created: %s" % (", ".join(missing_files))
            )

        committed_files: list[str] = []
        try:
            for output_file in self._output_files:
                if output_file.temporary:
                    fileutils.try_remove(os.path.join(self._temp, output_file.path))
                else:
                    destination = output_file.path
                    source = fileutils.reroot_path(self._temp, destination)

                    fileutils.move_file(source, destination)
                    committed_files.append(destination)
        except Exception:
            # Cleanup after failed commit;
            for fpath in committed_files:
                fileutils.try_remove(fpath)
            raise

        self._proc = None
        self._temp = None

    def to_call(self, temp: PathTypes) -> list[str]:
        return [self._to_path(temp, value) for value in self._command]

    @property
    def temp_dir(self) -> str | None:
        return self._temp

    @property
    def set_cwd(self) -> bool:
        return self._set_cwd

    @property
    def stdin(self) -> WrappedPipeType:
        return self._stdin

    @property
    def stdout(self) -> WrappedPipeType:
        return self._stdout

    @property
    def stderr(self) -> WrappedPipeType:
        return self._stderr

    @property
    def pid(self) -> int | None:
        return None if self._proc is None else self._proc.pid

    @property
    def terminated(self) -> bool:
        return self._terminated

    def _to_path(self, temp: PathTypes, value: str | AtomicFileTypes) -> str:
        if self._set_cwd:
            return self._to_abs_path(value)
        elif isinstance(value, Executable):
            executable = value.path
            if executable == "%(PYTHON)s":
                executable = sys.executable

            return executable
        elif isinstance(value, InputFile):
            if value.temporary:
                return os.path.join(temp, value.path)

            return value.path
        elif isinstance(value, OutputFile):
            return fileutils.reroot_path(temp, value.path)
        elif isinstance(value, AuxiliaryFile):
            return value.path
        else:
            return value.replace("%(TEMP_DIR)s", str(temp))

    def _to_abs_path(self, value: str | AtomicFileTypes) -> str:
        if isinstance(value, Executable):
            executable = value.path
            if executable == "%(PYTHON)s":
                executable = sys.executable
            elif os.path.dirname(executable) and os.path.relpath(executable):
                executable = os.path.abspath(executable)

            return executable
        elif isinstance(value, InputFile):
            if value.temporary:
                return value.path

            return os.path.abspath(value.path)
        elif isinstance(value, OutputFile):
            return os.path.basename(value.path)
        elif isinstance(value, AuxiliaryFile):
            return os.path.abspath(value.path)
        else:
            return value.replace("%(TEMP_DIR)s", ".")

    def _record_atomic_file(self, value: _AtomicFile) -> None:
        if isinstance(value, AuxiliaryFile):
            self._auxiliary_files.add(value.path)
        elif isinstance(value, Executable):
            self._executables.add(value.path)
        elif isinstance(value, InputFile):
            self._input_files.add(value)
        elif isinstance(value, OutputFile):
            basename = value.basename()
            # All output is saved in the same folder, so it is not possible to have
            # different output files with the same basename. It is however allowed to
            # use the same instance multiple times
            if basename in self._output_basenames and value not in self._output_files:
                raise CmdError(f"multiple output files with name {basename!r}")

            self._output_files.add(value)
            self._output_basenames.add(basename)
        else:
            raise TypeError(value)

    def _wrap_pipe(
        self,
        filetype: type[InputFile | OutputFile],
        pipe: PipeType,
        out_name: str | None = None,
    ) -> WrappedPipeType:
        if pipe is None:
            return self._wrap_pipe_default(filetype=filetype, out_name=out_name)
        elif isinstance(pipe, int):
            if pipe == AtomicCmd.DEVNULL or (
                pipe == AtomicCmd.PIPE and out_name == "stdout"
            ):
                return pipe
        elif isinstance(pipe, _AtomicFile):
            if isinstance(pipe, filetype):
                return pipe

            raise ValueError(f"expected {filetype.__name__}, but got {pipe}")
        elif isinstance(pipe, AtomicCmd):
            # Piping with a AtomicCmd is only allowed for STDIN
            if out_name is None:
                return pipe
        elif isinstance(pipe, (str, Path)):
            return filetype(pipe)

        raise ValueError(pipe)

    def _wrap_pipe_default(
        self,
        filetype: type[InputFile | OutputFile],
        out_name: str | None = None,
    ) -> WrappedPipeType:
        # Pipes without a default destination (stdin) are explicitly closed
        if out_name is None:
            return self.DEVNULL

        assert isinstance(self._command[0], Executable)
        executable = self._command[0].path
        if executable == "%(PYTHON)s":
            executable = sys.executable

        executable = os.path.basename(executable)

        # TODO: Use more sensible filename, e.g. log_{exec}_{counter}.{stdout/err}
        filename = "pipe_%s_%i.%s" % (executable, id(self), out_name)

        return filetype(filename, temporary=True)

    @classmethod
    def _open_pipe(
        cls,
        temp_dir: PathTypes,
        pipe: WrappedPipeType,
        mode: str,
    ) -> int | IO[bytes] | None:
        if isinstance(pipe, int):
            return pipe
        elif isinstance(pipe, AtomicCmd):
            if pipe._proc is None:
                raise CmdError("attempted to pipe non-running command")
            elif pipe._proc.stdout is None:
                raise CmdError("attempted to pipe from command without stdout=PIPE")

            return pipe._proc.stdout

        if pipe.temporary or isinstance(pipe, OutputFile):
            return open(fileutils.reroot_path(temp_dir, pipe.path), mode)  # noqa: SIM115

        return open(pipe.path, mode)  # noqa: SIM115

    def __enter__(self) -> Self:
        return self

    def __exit__(self, typ: object, exc: object, tb: object) -> None:
        self.terminate()
        self.join()

    def __str__(self) -> str:
        return pformat(self)


class _CommandSet:
    def __init__(self, commands: Iterable[AtomicCmd | _CommandSet]) -> None:
        self._commands: tuple[CommandTypes, ...] = safe_coerce_to_tuple(commands)
        if not self._commands:
            raise CmdError("Empty list passed to command set")

        self._validate_commands()

    def commit(self) -> None:
        committed_files: set[str] = set()
        try:
            for command in self._commands:
                command.commit()
                committed_files.update(command.output_files)
        except Exception:
            # Cleanup after failed commit; files belonging to the command currently
            # being comitted are cleaned up in `AtomicCmd.commit`
            for fpath in committed_files:
                fileutils.try_remove(fpath)
            raise

    @property
    def input_files(self) -> set[str]:
        return {it for cmd in self._commands for it in cmd.input_files}

    @property
    def output_files(self) -> set[str]:
        return {it for cmd in self._commands for it in cmd.output_files}

    @property
    def auxiliary_files(self) -> set[str]:
        return {it for cmd in self._commands for it in cmd.auxiliary_files}

    @property
    def executables(self) -> set[str]:
        return {it for cmd in self._commands for it in cmd.executables}

    @property
    def requirements(self) -> set[Requirement]:
        return {it for cmd in self._commands for it in cmd.requirements}

    @property
    def expected_temp_files(self) -> set[str]:
        return {it for cmd in self._commands for it in cmd.expected_temp_files}

    @property
    def optional_temp_files(self) -> set[str]:
        return {it for cmd in self._commands for it in cmd.optional_temp_files}

    def terminate(self) -> None:
        for command in self._commands:
            command.terminate()

    def _validate_commands(self) -> None:
        if len(self._commands) != len(set(self._commands)):
            raise ValueError(
                f"Same command included multiple times in {self.__class__.__name__}"
            )

        filenames: dict[str, int] = collections.defaultdict(int)
        for command in self._commands:
            for filename in command.expected_temp_files:
                filenames[filename] += 1
            for filename in command.optional_temp_files:
                filenames[filename] += 1

        clobbered = [filename for (filename, count) in filenames.items() if (count > 1)]
        if any(clobbered):
            raise CmdError(
                "Commands clobber each others' files: {}".format(", ".join(clobbered))
            )


class ParallelCmds(_CommandSet):
    """This class wraps a set of AtomicCmds, running them in parallel.
    This corresponds to a set of piped commands, which only terminate
    when all parts of the pipe have terminated. For example:
    $ dmesg | grep -i segfault | gzip > log.txt.gz

    In case of any one sub-command failing, the remaining commands are
    automatically terminated. This is done to ensure that commands waiting
    on pipes are not left running indefinitely.

    Note that only AtomicCmds are allowed as sub-commands for this class,
    since the model requires non-blocking commands."""

    def __init__(self, commands: Iterable[AtomicCmd]) -> None:
        self._joinable = False

        commands = tuple(commands)
        for command in commands:
            if not isinstance(command, AtomicCmd):
                raise CmdError("ParallelCmds must only contain AtomicCmds")
        _CommandSet.__init__(self, commands)

    def run(self, temp: PathTypes) -> None:
        for command in self._commands:
            command.run(temp)
        self._joinable = True

    def ready(self) -> bool:
        return all(cmd.ready() for cmd in self._commands)

    def join(self) -> JoinType:
        sleep_time = 0.05
        commands: list[tuple[int, AtomicCmd]] = []
        for idx, command in enumerate(self._commands):
            if not isinstance(command, AtomicCmd):
                raise TypeError("non AtomicCmd in ParallelCmds")

            commands.append((idx, command))

        return_codes: JoinType = [None] * len(commands)
        if not self._joinable:
            return return_codes

        while commands and not any(return_codes):
            try:
                # Wait for arbitrary command
                commands[0][1].join(sleep_time if len(commands) > 1 else None)
            except subprocess.TimeoutExpired:
                sleep_time = min(1, sleep_time * 2)

            for index, command in tuple(commands):
                if command.ready():
                    (return_code,) = command.join()
                    return_codes[index] = return_code
                    commands.remove((index, command))
                    sleep_time = 0.05

        if any(return_codes):
            for index, command in commands:
                command.terminate()
                (return_codes[index],) = command.join()

        return return_codes

    def __str__(self) -> str:
        return pformat(self)


class SequentialCmds(_CommandSet):
    """This class wraps a set of AtomicCmds, running them sequentially.
    This class therefore corresponds a set of lines in a bash script,
    each of which invokes a foreground job. For example:
    $ bcftools view snps.bcf | bgzip > snps.vcf.bgz
    $ tabix snps.vcf.bgz

    The list of commands may include any type of command. Note that
    the run function only returns once each sub-command has completed.
    A command is only executed if the previous command in the sequence
    was successfully completed, and as a consequence the return codes
    of a failed SequentialCommand may contain None."""

    def __init__(self, commands: Iterable[AtomicCmd | _CommandSet]) -> None:
        self._ready = False

        commands = safe_coerce_to_tuple(commands)
        for command in commands:
            if not isinstance(command, (AtomicCmd, _CommandSet)):
                raise CmdError(
                    "ParallelCmds must only contain AtomicCmds or other ParallelCmds!"
                )
        _CommandSet.__init__(self, commands)

    def run(self, temp: PathTypes) -> None:
        self._ready = False
        for command in self._commands:
            command.run(temp)
            if any(command.join()):
                break

        self._ready = True

    def ready(self) -> bool:
        return self._ready

    def join(self) -> JoinType:
        return_codes: JoinType = []
        for command in self._commands:
            return_codes.extend(command.join())

        return return_codes

    def __str__(self) -> str:
        return pformat(self)


CommandTypes = Union[AtomicCmd, ParallelCmds, SequentialCmds]


def _describe_cls(command: ParallelCmds | SequentialCmds) -> str:
    if isinstance(command, ParallelCmds):
        return "Parallel processes"

    return "Sequential processes"


def _collect_stats(
    command: CommandTypes,
    ids: dict[CommandTypes, int],
    pipes: dict[AtomicCmd, AtomicCmd],
) -> None:
    assert command not in ids

    if isinstance(command, AtomicCmd):
        ids[command] = len(ids) + 1
        if isinstance(command.stdin, AtomicCmd):
            pipes[command.stdin] = command
    elif isinstance(command, (ParallelCmds, SequentialCmds)):
        for subcmd in command._commands:
            _collect_stats(subcmd, ids, pipes)
    else:
        raise ValueError(command)


def _build_status(command: AtomicCmd, indent: int, lines: list[str]) -> None:
    prefix = " " * indent + "Status  = "
    if command._proc:
        if command.ready():
            return_code = tuple(command.join())
            if command.terminated:
                lines.append(prefix + "Automatically terminated by PALEOMIX")
            elif isinstance(return_code[0], str):
                lines.append(prefix + "Terminated with signal %s" % return_code)
            else:
                lines.append(prefix + "Exited with return-code %i" % return_code)
        else:
            lines.append(prefix + "Running")


def _build_stdin(
    command: AtomicCmd,
    ids: dict[CommandTypes, int],
    indent: int,
    lines: list[str],
) -> None:
    pipe = command.stdin
    prefix = "{}{}   = ".format(" " * indent, "STDIN")
    if isinstance(pipe, AtomicCmd):
        if pipe in ids:
            lines.append("%sPiped from process %i" % (prefix, ids[pipe]))
        else:
            lines.append(f"{prefix}<PIPE>")
    elif isinstance(pipe, (InputFile, TempInputFile)):
        temp = "${TEMP_DIR}" if command.temp_dir is None else command.temp_dir
        path = command._to_path(temp, pipe)
        lines.append(f"{prefix}{shlex.quote(path)}")


def _build_stdout(
    command: AtomicCmd,
    ids: dict[CommandTypes, int],
    pipes: dict[AtomicCmd, AtomicCmd],
    indent: int,
    lines: list[str],
) -> None:
    prefix = "{}STDOUT  = ".format(" " * indent)

    pipe = command.stdout
    if command in pipes:
        pipe = pipes[command]
        lines.append("%sPiped to process %i" % (prefix, ids[pipe]))
    elif isinstance(pipe, (OutputFile, TempOutputFile)):
        temp = "${TEMP_DIR}" if command.temp_dir is None else command.temp_dir
        path = command._to_path(temp, pipe)

        lines.append(f"{prefix}{shlex.quote(path)}")
    elif pipe == command.PIPE:
        lines.append(f"{prefix}<PIPE>")
    elif pipe == command.DEVNULL:
        lines.append(f"{prefix}/dev/null")
    else:
        raise TypeError(pipe)  # pragma: no coverage


def _build_stderr(command: AtomicCmd, indent: int, lines: list[str]) -> None:
    prefix = "{}STDERR  = ".format(" " * indent)

    pipe = command.stderr
    if isinstance(pipe, (OutputFile, TempOutputFile)):
        temp = "${TEMP_DIR}" if command.temp_dir is None else command.temp_dir
        path = command._to_path(temp, pipe)

        lines.append(f"{prefix}{shlex.quote(path)}")
    elif pipe == command.DEVNULL:
        lines.append(f"{prefix}/dev/null")
    else:
        raise TypeError(pipe)  # pragma: no coverage


def _build_cwd(command: AtomicCmd, indent: int, lines: list[str]) -> None:
    prefix = " " * indent + "CWD     = "
    if command.temp_dir:
        if command.set_cwd:
            lines.append(f"{prefix}{shlex.quote(command.temp_dir)}")
        else:
            lines.append(f"{prefix}{shlex.quote(os.getcwd())}")
    elif command.set_cwd:
        lines.append("{}{}".format(prefix, shlex.quote("${TEMP_DIR}")))


def _pformat(
    command: CommandTypes,
    *,
    ids: dict[CommandTypes, int],
    pipes: dict[AtomicCmd, AtomicCmd],
    indent: int,
    lines: list[str],
    include_prefix: bool = True,
) -> None:
    s_prefix = ""
    if include_prefix:
        s_prefix = " " * indent
        if isinstance(command, AtomicCmd):
            cmd_id = ids[command]
            lines.append(s_prefix + "Process %i:" % (cmd_id,))
            s_prefix += "  "
    s_prefix_len = len(s_prefix)

    if isinstance(command, AtomicCmd):
        temp = command.temp_dir or "${TEMP_DIR}"

        c_prefix = s_prefix + "Command = "
        for line in _pformat_list(command.to_call(temp)).split("\n"):
            lines.append(f"{c_prefix}{line}")
            c_prefix = " " * len(c_prefix)

        _build_status(command, s_prefix_len, lines)
        _build_stdin(command, ids, s_prefix_len, lines)
        _build_stdout(command, ids, pipes, s_prefix_len, lines)
        _build_stderr(command, s_prefix_len, lines)
        _build_cwd(command, s_prefix_len, lines)
    else:
        lines.append(f"{s_prefix}{_describe_cls(command)}:")
        for subcmd_idx, subcmd in enumerate(command._commands):
            if subcmd_idx:
                lines.append("")

            _pformat(subcmd, ids=ids, pipes=pipes, indent=s_prefix_len + 2, lines=lines)


def _pformat_list(lst: list[str], width: int = 80) -> str:
    """Return a printable representation of a list, where line-breaks
    are inserted between items to minimize the number of lines with a
    width greater than 'width'. Very long items may cause this maximum
    to be exceeded."""
    result: list[list[str]] = [[]]
    current_width = 0
    for item in (shlex.quote(value) for value in lst):
        if current_width + len(item) + 1 > width:
            if not result[-1]:
                result[-1] = [item]
            else:
                result.append([item])

            current_width = len(item) + 1
        else:
            result[-1].append(item)
            current_width += len(item) + 1

    return " \\\n    ".join(" ".join(line) for line in result)


def pformat(command: CommandTypes) -> str:
    """Returns a human readable description of an Atomic Cmd or Atomic Set
    of commands. This is currently equivalent to str(cmd_obj)."""
    if not isinstance(command, (AtomicCmd, ParallelCmds, SequentialCmds)):
        raise TypeError(command)

    lines: list[str] = []
    ids: dict[CommandTypes, int] = {}
    pipes: dict[AtomicCmd, AtomicCmd] = {}

    _collect_stats(command, ids, pipes)
    _pformat(command, ids=ids, pipes=pipes, indent=0, lines=lines, include_prefix=False)

    return "\n".join(lines)
