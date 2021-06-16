#!/usr/bin/env python3
from __future__ import annotations

import atexit
import collections
import os
import shlex
import signal
import subprocess
import sys
import time
import weakref
from pathlib import Path
from typing import (
    IO,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    NoReturn,
    Optional,
    Set,
    Tuple,
    Type,
    Union,
)

import paleomix.common.fileutils as fileutils
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.versions import Requirement


class CmdError(RuntimeError):
    """Exception raised for AtomicCmd specific errors."""

    def __init__(self, msg: object):
        RuntimeError.__init__(self, msg)


class _AtomicFile:
    def __init__(self, path: Union[Path, str]):
        self.path = fileutils.fspath(path)

        if isinstance(self.path, bytes):
            raise ValueError(path)

    def basename(self):
        return os.path.basename(self.path)

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.path)


class AuxilleryFile(_AtomicFile):
    pass


class Executable(_AtomicFile):
    pass


class _IOFile(_AtomicFile):
    def __init__(self, path: Union[Path, str], temporary: bool = False):
        super().__init__(path)
        self.temporary = bool(temporary)

        if temporary and os.path.dirname(self.path):
            raise ValueError(self.path)

    def __repr__(self):
        return "%s(%r, %s)" % (self.__class__.__name__, self.path, self.temporary)


class InputFile(_IOFile):
    pass


class TempInputFile(InputFile):
    def __init__(self, path: Union[Path, str]):
        super().__init__(os.path.basename(path), temporary=True)


class OutputFile(_IOFile):
    pass


class TempOutputFile(OutputFile):
    def __init__(self, path: str):
        super().__init__(os.path.basename(path), temporary=True)


# Possible types of .stdin/.stdout/.stderr, int being either DEVNULL or PIPE
WrappedPipeType = Union[int, _IOFile, "AtomicCmd"]
# Types that can be passed as values for STDIN, STDOUT, and STDERR
PipeType = Union[None, str, Path, WrappedPipeType]

# Pos
OptionValueType = Union[str, float, _IOFile, None]
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

    InFile = InputFile
    OutFile = OutputFile
    AuxFile = AuxilleryFile
    Executable = Executable

    _command: List[Union[str, _AtomicFile]]
    # _proc:  = None
    _temp: Optional[str]
    _running: bool
    _set_cwd: bool
    _terminated: bool

    _executables: Set[str]
    _requirements: Set[Requirement]
    _input_files: Set[_IOFile]
    _output_files: Set[_IOFile]
    _output_basenames: Set[str]
    _auxiliary_files: Set[str]

    _stdin: WrappedPipeType
    _stdout: WrappedPipeType
    _stderr: WrappedPipeType

    def __init__(
        self,
        command: Any,
        stdin: PipeType = None,
        stdout: PipeType = None,
        stderr: PipeType = None,
        set_cwd: bool = False,
        extra_files: Iterable[_AtomicFile] = (),
        requirements: Iterable[Requirement] = (),
    ):
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
            raise ValueError(self._command[0])

        self._stdin = stdin = self._wrap_pipe(InputFile, stdin)
        self._stdout = stdout = self._wrap_pipe(OutputFile, stdout, "stdout")
        self._stderr = stderr = self._wrap_pipe(OutputFile, stderr, "stderr")

        pipes = []  # type: List[_AtomicFile]
        for pipe in (stdin, stdout, stderr):
            if isinstance(pipe, _AtomicFile):
                pipes.append(pipe)

        self.add_extra_files(extra_files)
        self.add_extra_files(pipes)

        for value in self._requirements:
            if not isinstance(value, Requirement):
                raise TypeError(value)

    def append(self, *args: Any):
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

    def add_extra_files(self, files: Iterable[Any]):
        if self._proc is not None:
            raise CmdError("cannot modify already started command")

        for value in files:
            if not isinstance(value, _AtomicFile):
                raise ValueError(value)

            self._record_atomic_file(value)

    def append_options(
        self,
        options: OptionsType,
        pred: Callable[[str], bool] = lambda s: s.startswith("-"),
    ):
        if not isinstance(options, dict):
            raise TypeError("options must be dict, not {!r}".format(options))

        for (key, values) in options.items():
            if not isinstance(key, str):
                raise ValueError("keys must be strings, not %r" % (key,))
            elif not pred(key):
                continue

            if isinstance(values, (list, tuple)):
                for value in values:
                    if not isinstance(value, (int, str, float, _AtomicFile)):
                        raise ValueError(value)

                    self.append(key)
                    self.append(value)
            elif values is None:
                self.append(key)
            elif isinstance(values, (int, str, float, _AtomicFile)):
                self.append(key)
                self.append(values)
            else:
                raise ValueError(values)

    def merge_options(
        self,
        user_options: OptionsType,
        fixed_options: OptionsType = {},
        blacklisted_options: Iterable[str] = (),
        pred: Callable[[str], bool] = lambda s: s.startswith("-"),
    ):
        if not isinstance(fixed_options, dict):
            raise TypeError("options must be dict, not {!r}".format(fixed_options))
        elif not isinstance(user_options, dict):
            raise TypeError("user_options must be dict, not {!r}".format(user_options))

        errors: List[str] = []
        for key in user_options.keys() & fixed_options.keys():
            errors.append("{} cannot be overridden".format(key))

        for key in user_options.keys() & blacklisted_options:
            errors.append("{} is not supported".format(key))

        if errors:
            raise CmdError(
                "invalid command-line options for {!r}: {}".format(
                    " ".join(self.to_call("%(TEMP_DIR)s")), "\n".join(errors)
                )
            )

        self.append_options(fixed_options, pred=pred)
        self.append_options(user_options, pred=pred)

    @property
    def output_files(self):
        return frozenset(
            afile.path for afile in self._output_files if not afile.temporary
        )

    @property
    def expected_temp_files(self):
        return frozenset(
            afile.basename() for afile in self._output_files if not afile.temporary
        )

    @property
    def optional_temp_files(self):
        return frozenset(
            afile.basename() for afile in self._output_files if afile.temporary
        )

    @property
    def input_files(self):
        return frozenset(
            afile.path for afile in self._input_files if not afile.temporary
        )

    @property
    def executables(self):
        return frozenset(self._executables)

    @property
    def auxiliary_files(self):
        return frozenset(self._auxiliary_files)

    @property
    def requirements(self):
        return self._requirements

    def run(self, temp: Union[Path, str]):
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

            self._proc = subprocess.Popen(
                call,
                stdin=stdin,
                stdout=stdout,
                stderr=stderr,
                cwd=cwd,
                start_new_session=True,
            )
        except Exception as error:
            message = "Error running commands:\n  Call = %r\n  Error = %r"
            raise CmdError(message % (self._command, error))
        finally:
            # Close pipes to allow the command to recieve SIGPIPE
            for handle in (stdin, stdout, stderr):
                if not (handle is None or isinstance(handle, int)):
                    handle.close()

        # Allow subprocesses to be killed in case of a SIGTERM
        _add_to_killlist(self._proc)

    def ready(self):
        """Returns true if the command has been run to completion,
        regardless of wether or not an error occured."""
        return self._proc and self._proc.poll() is not None

    def join(self) -> JoinType:
        """Similar to Popen.wait(), but returns the value wrapped in a list.
        Must be called before calling commit."""
        if not self._proc:
            return [None]

        return_code = self._proc.wait()
        self._running = False
        if return_code < 0:
            return_code = signal.Signals(-return_code).name
        return [return_code]

    def wait(self):
        """Equivalent to Subproces.wait. This function should only
        be used in contexts where a AtomicCmd needs to be combined
        with Subprocesses, as it does not exist for AtomicSets."""
        return self.join()[0]

    def terminate(self):
        """Sends SIGTERM to process if it is still running.
        Has no effect if the command has already finished."""
        if self._proc and self._proc.poll() is None:
            try:
                os.killpg(self._proc.pid, signal.SIGTERM)
                self._terminated = True
            except OSError:
                pass  # Already dead / finished process

    def commit(self, temp: Union[Path, str]):
        temp = fileutils.fspath(temp)
        if not self.ready():
            raise CmdError("Attempting to commit before command has completed")

        assert self._temp is not None
        if self._running:
            raise CmdError("Called 'commit' before calling 'join'")
        elif not os.path.samefile(self._temp, temp):
            raise CmdError(
                "Mismatch between previous and current temp folders"
                ": %r != %s" % (self._temp, temp)
            )

        expected_files = set(os.path.basename(fpath) for fpath in self.output_files)
        missing_files = expected_files - set(os.listdir(temp))
        if missing_files:
            raise CmdError(
                "Expected files not created: %s" % (", ".join(missing_files))
            )

        committed_files: List[str] = []
        try:
            for output_file in self._output_files:
                if output_file.temporary:
                    fileutils.try_remove(os.path.join(temp, output_file.path))
                else:
                    destination = output_file.path
                    source = fileutils.reroot_path(temp, destination)

                    fileutils.move_file(source, destination)
                    committed_files.append(destination)
        except Exception:
            # Cleanup after failed commit
            for fpath in committed_files:
                fileutils.try_remove(fpath)
            raise

        self._proc = None
        self._temp = None

    def to_call(self, temp: str) -> List[str]:
        return [self._to_path(temp, value) for value in self._command]

    def _to_path(self, temp: str, value: Any) -> str:
        if self._set_cwd:
            if isinstance(value, Executable):
                executable = value.path
                if executable == "%(PYTHON)s":
                    executable = sys.executable
                elif os.path.dirname(executable) and os.path.relpath(executable):
                    return os.path.abspath(executable)

                return executable
            elif isinstance(value, InputFile):
                if value.temporary:
                    return value.path

                return os.path.abspath(value.path)
            elif isinstance(value, OutputFile):
                return os.path.basename(value.path)
            elif isinstance(value, AuxilleryFile):
                return os.path.abspath(value.path)
            else:
                return value.replace("%(TEMP_DIR)s", ".")
        else:
            if isinstance(value, Executable):
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
            elif isinstance(value, AuxilleryFile):
                return value.path
            else:
                return value.replace("%(TEMP_DIR)s", temp)

    def _record_atomic_file(self, value: _AtomicFile):
        if isinstance(value, AuxilleryFile):
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
                raise CmdError("multiple output files with name %r" % (basename,))

            self._output_files.add(value)
            self._output_basenames.add(basename)
        else:
            raise ValueError(value)

    def _wrap_pipe(
        self,
        filetype: Union[Type[InputFile], Type[OutputFile]],
        pipe: PipeType,
        out_name: Optional[str] = None,
    ) -> WrappedPipeType:
        if pipe is None:
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
        elif isinstance(pipe, int):
            if pipe == AtomicCmd.DEVNULL or (
                pipe == AtomicCmd.PIPE and out_name == "stdout"
            ):
                return pipe
        elif isinstance(pipe, _AtomicFile):
            if isinstance(pipe, filetype):
                return pipe

            raise ValueError("expected %s, but got %s" % (filetype, pipe))
        elif isinstance(pipe, AtomicCmd):
            # Piping with a AtomicCmd is only allowed for STDIN
            if out_name is None:
                return pipe
        elif isinstance(pipe, (str, Path)):
            return filetype(pipe)

        raise ValueError(pipe)

    @classmethod
    def _open_pipe(
        cls,
        temp_dir: str,
        pipe: WrappedPipeType,
        mode: str,
    ) -> Union[int, IO[bytes], None]:
        if isinstance(pipe, int):
            return pipe
        elif isinstance(pipe, AtomicCmd):
            if pipe._proc is None:
                raise CmdError("attempted to pipe non-running command")
            elif pipe._proc.stdout is None:
                raise CmdError("attempted to pipe from command without stdout=PIPE")

            return pipe._proc.stdout

        if pipe.temporary or isinstance(pipe, OutputFile):
            return open(fileutils.reroot_path(temp_dir, pipe.path), mode)

        return open(pipe.path, mode)

    def __enter__(self):
        return self

    def __exit__(self, type: Any, _value: Any, _traceback: Any):
        self.terminate()
        self.join()

    def __str__(self) -> str:
        return pformat(self)


class _CommandSet:
    def __init__(self, commands: Iterable[Union[AtomicCmd, "_CommandSet"]]) -> None:
        self._commands: Tuple[CommandTypes, ...] = safe_coerce_to_tuple(commands)
        if not self._commands:
            raise CmdError("Empty list passed to command set")

        self._validate_commands()

    def commit(self, temp: str) -> None:
        committed_files: Set[str] = set()
        try:
            for command in self._commands:
                command.commit(temp)
                committed_files.update(command.output_files)
        except Exception:
            # Cleanup after failed commit
            for fpath in committed_files:
                fileutils.try_remove(fpath)
            raise

    @property
    def input_files(self) -> Set[str]:
        return set(it for cmd in self._commands for it in cmd.input_files)

    @property
    def output_files(self) -> Set[str]:
        return set(it for cmd in self._commands for it in cmd.output_files)

    @property
    def auxiliary_files(self) -> Set[str]:
        return set(it for cmd in self._commands for it in cmd.auxiliary_files)

    @property
    def executables(self) -> Set[str]:
        return set(it for cmd in self._commands for it in cmd.executables)

    @property
    def requirements(self) -> Set[Requirement]:
        return set(it for cmd in self._commands for it in cmd.requirements)

    @property
    def expected_temp_files(self) -> Set[str]:
        return set(it for cmd in self._commands for it in cmd.expected_temp_files)

    @property
    def optional_temp_files(self) -> Set[str]:
        return set(it for cmd in self._commands for it in cmd.optional_temp_files)

    @property
    def stdout(self):
        raise CmdError(
            "%s does not implement property 'stdout'!" % (self.__class__.__name__,)
        )

    def terminate(self) -> None:
        for command in self._commands:
            command.terminate()

    def _validate_commands(self):
        if len(self._commands) != len(set(self._commands)):
            raise ValueError(
                "Same command included multiple times in %s"
                % (self.__class__.__name__,)
            )

        filenames: Dict[str, int] = collections.defaultdict(int)
        for command in self._commands:
            for filename in command.expected_temp_files:
                filenames[filename] += 1
            for filename in command.optional_temp_files:
                filenames[filename] += 1

        clobbered = [filename for (filename, count) in filenames.items() if (count > 1)]
        if any(clobbered):
            raise CmdError(
                "Commands clobber each others' files: %s" % (", ".join(clobbered),)
            )


class ParallelCmds(_CommandSet):
    """This class wraps a set of AtomicCmds, running them in parallel.
    This corresponds to a set of piped commands, which only terminate
    when all parts of the pipe have terminated. For example:
    $ dmesg | grep -i segfault | gzip > log.txt.gz

    In case of any one sub-command failing, the remaining commands are
    automatically terminated. This is done to ensure that commands waiting
    on pipes are not left running indefinetly.

    Note that only AtomicCmds and ParallelCmds are allowed as
    sub-commands for this class, since the model requires non-
    blocking commands."""

    def __init__(self, commands: Iterable[Union[AtomicCmd, "ParallelCmds"]]):
        self._joinable = False

        commands = tuple(commands)
        for command in commands:
            if not isinstance(command, (AtomicCmd, ParallelCmds)):
                raise CmdError(
                    "ParallelCmds must only contain AtomicCmds or other ParallelCmds!"
                )
        _CommandSet.__init__(self, commands)

    def run(self, temp: str) -> None:
        for command in self._commands:
            command.run(temp)
        self._joinable = True

    def ready(self) -> bool:
        return all(cmd.ready() for cmd in self._commands)

    def join(self) -> JoinType:
        sleep_time = 0.05
        commands = list(enumerate(self._commands))
        return_codes: List[JoinType] = [[None]] * len(commands)
        while commands and self._joinable:
            for (index, command) in list(commands):
                if command.ready():
                    return_codes[index] = command.join()
                    commands.remove((index, command))
                    sleep_time = 0.05
                elif any(any(codes) for codes in return_codes):
                    command.terminate()
                    return_codes[index] = command.join()
                    commands.remove((index, command))
                    sleep_time = 0.05

            time.sleep(sleep_time)
            sleep_time = min(1, sleep_time * 2)

        result: JoinType = []
        for return_code in return_codes:
            result.extend(return_code)

        return result

    def __str__(self) -> str:
        return pformat(self)


class SequentialCmds(_CommandSet):
    """This class wraps a set of AtomicCmds, running them sequentially.
    This class therefore corresponds a set of lines in a bash script,
    each of which invokes a forground job. For example:
    $ bcftools view snps.bcf | bgzip > snps.vcf.bgz
    $ tabix snps.vcf.bgz

    The list of commands may include any type of command. Note that
    the run function only returns once each sub-command has completed.
    A command is only executed if the previous command in the sequence
    was succesfully completed, and as a consequence the return codes
    of a failed SequentialCommand may contain None."""

    def __init__(self, commands: Iterable[Union[AtomicCmd, _CommandSet]]):
        self._ready = False

        commands = safe_coerce_to_tuple(commands)
        for command in commands:
            if not isinstance(command, (AtomicCmd, _CommandSet)):
                raise CmdError(
                    "ParallelCmds must only contain AtomicCmds or other ParallelCmds!"
                )
        _CommandSet.__init__(self, commands)

    def run(self, temp: str) -> None:
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


def _describe_cls(command: Union[ParallelCmds, SequentialCmds]) -> str:
    if isinstance(command, ParallelCmds):
        return "Parallel processes"
    elif isinstance(command, SequentialCmds):
        return "Sequential processes"
    assert False  # pragma: no coverage


def _collect_stats(
    command: CommandTypes,
    ids: Dict[CommandTypes, int],
    pipes: Dict[AtomicCmd, AtomicCmd],
):
    assert command not in ids

    if isinstance(command, AtomicCmd):
        ids[command] = len(ids) + 1
        if isinstance(command._stdin, AtomicCmd):
            pipes[command._stdin] = command
    elif isinstance(command, (ParallelCmds, SequentialCmds)):
        for subcmd in command._commands:
            _collect_stats(subcmd, ids, pipes)
    else:
        assert False  # pragma: no coverage


def _build_status(command: AtomicCmd, indent: int, lines: List[str]) -> None:
    prefix = " " * indent + "Status  = "
    if command._proc:
        if command.ready():
            return_code = tuple(command.join())
            if command._terminated:
                lines.append(prefix + "Automatically terminated by PALEOMIX")
            elif isinstance(return_code[0], str):
                lines.append(prefix + "Terminated with signal %s" % return_code)
            else:
                lines.append(prefix + "Exited with return-code %i" % return_code)
        else:
            lines.append(prefix + "Running")


def _build_stdin(
    command: AtomicCmd,
    ids: Dict[CommandTypes, int],
    indent: int,
    lines: List[str],
) -> None:
    pipe = command._stdin
    prefix = "%s%s   = " % (" " * indent, "STDIN")
    if isinstance(pipe, AtomicCmd):
        if pipe in ids:
            lines.append("%sPiped from process %i" % (prefix, ids[pipe]))
        else:
            lines.append("%s<PIPE>" % (prefix,))
    elif isinstance(pipe, (InputFile, TempInputFile)):
        temp = "${TEMP_DIR}" if command._temp is None else command._temp
        path = command._to_path(temp, pipe)
        lines.append("%s%s" % (prefix, shlex.quote(path)))


def _build_stdout(
    command: AtomicCmd,
    ids: Dict[CommandTypes, int],
    pipes: Dict[AtomicCmd, AtomicCmd],
    indent: int,
    lines: List[str],
) -> None:
    prefix = "%sSTDOUT  = " % (" " * indent,)

    pipe = command._stdout
    if command in pipes:
        pipe = pipes[command]
        lines.append("%sPiped to process %i" % (prefix, ids[pipe]))
    elif isinstance(pipe, (OutputFile, TempOutputFile)):
        temp = "${TEMP_DIR}" if command._temp is None else command._temp
        path = command._to_path(temp, pipe)

        lines.append("%s%s" % (prefix, shlex.quote(path)))
    elif pipe == command.PIPE:
        lines.append("%s<PIPE>" % (prefix,))
    elif pipe == command.DEVNULL:
        lines.append("%s/dev/null" % (prefix,))
    else:
        assert False, pipe  # pragma: no coverage


def _build_stderr(command: AtomicCmd, indent: int, lines: List[str]) -> None:
    prefix = "%sSTDERR  = " % (" " * indent,)

    pipe = command._stderr
    if isinstance(pipe, (OutputFile, TempOutputFile)):
        pipe = pipe
        temp = "${TEMP_DIR}" if command._temp is None else command._temp
        path = command._to_path(temp, pipe)

        lines.append("%s%s" % (prefix, shlex.quote(path)))
    elif pipe == command.DEVNULL:
        lines.append("%s/dev/null" % (prefix,))
    else:
        assert False, pipe  # pragma: no coverage


def _build_cwd(command: AtomicCmd, indent: int, lines: List[str]) -> None:
    prefix = " " * indent + "CWD     = "
    if command._temp:
        if command._set_cwd:
            lines.append("%s%s" % (prefix, shlex.quote(command._temp)))
        else:
            lines.append("%s%s" % (prefix, shlex.quote(os.getcwd())))
    elif command._set_cwd:
        lines.append("%s%s" % (prefix, shlex.quote("${TEMP_DIR}")))


def _pformat(
    command: CommandTypes,
    ids: Dict[CommandTypes, int],
    pipes: Dict[AtomicCmd, AtomicCmd],
    indent: int,
    lines: List[str],
    include_prefix: bool = True,
):
    s_prefix = ""
    if include_prefix:
        s_prefix = " " * indent
        if isinstance(command, AtomicCmd):
            cmd_id = ids[command]
            lines.append(s_prefix + "Process %i:" % (cmd_id,))
            s_prefix += "  "
    s_prefix_len = len(s_prefix)

    if isinstance(command, AtomicCmd):
        temp = command._temp or "${TEMP_DIR}"

        c_prefix = s_prefix + "Command = "
        for line in _pformat_list(command.to_call(temp)).split("\n"):
            lines.append("%s%s" % (c_prefix, line))
            c_prefix = " " * len(c_prefix)

        _build_status(command, s_prefix_len, lines)
        _build_stdin(command, ids, s_prefix_len, lines)
        _build_stdout(command, ids, pipes, s_prefix_len, lines)
        _build_stderr(command, s_prefix_len, lines)
        _build_cwd(command, s_prefix_len, lines)
    else:
        lines.append("%s%s:" % (s_prefix, _describe_cls(command)))
        for subcmd_idx, subcmd in enumerate(command._commands):
            if subcmd_idx:
                lines.append("")

            _pformat(subcmd, ids, pipes, s_prefix_len + 2, lines)


def _pformat_list(lst: List[Any], width: int = 80):
    """Return a printable representation of a list, where line-breaks
    are inserted between items to minimize the number of lines with a
    width greater than 'width'. Very long items may cause this maximum
    to be exceeded."""
    result: List[List[str]] = [[]]
    current_width = 0
    for item in (shlex.quote(str(value)) for value in lst):
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

    lines: List[str] = []
    ids: Dict[CommandTypes, int] = {}
    pipes: Dict[AtomicCmd, AtomicCmd] = {}

    _collect_stats(command, ids, pipes)
    _pformat(command, ids, pipes, 0, lines, False)

    return "\n".join(lines)


# The following ensures proper cleanup of child processes, for example in the
# case where multiprocessing.Pool.terminate() is called or if the script exits due to
# an unhandled exception.
_PROCS: List[weakref.ReferenceType[subprocess.Popen[Any]]] = []


@atexit.register
def _cleanup_children() -> None:
    for proc_ref in list(_PROCS):
        proc = proc_ref()
        try:
            if proc:
                proc.terminate()
        except OSError:
            # Ignore already closed processes, etc.
            pass


def _on_sig_term(signum: int, _frame: Any) -> NoReturn:
    _cleanup_children()
    sys.exit(-signum)


def _add_to_killlist(proc: subprocess.Popen[Any]) -> None:
    _PROCS.append(weakref.ref(proc, _PROCS.remove))


signal.signal(signal.SIGTERM, _on_sig_term)
