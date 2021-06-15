#!/usr/bin/env python3
import io
import os
import signal
import subprocess
import sys
import weakref
from typing import (
    IO,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Set,
    Tuple,
    Type,
    Union,
)

import paleomix.atomiccmd.pprint as atomicpp
import paleomix.common.fileutils as fileutils
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.versions import Requirement


class CmdError(RuntimeError):
    """Exception raised for AtomicCmd specific errors."""

    def __init__(self, msg: object):
        RuntimeError.__init__(self, msg)


class _AtomicFile:
    def __init__(self, path: str):
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
    def __init__(self, path: str, temporary: bool = False):
        super().__init__(path)
        self.temporary = bool(temporary)

        if temporary and os.path.dirname(self.path):
            raise ValueError(self.path)

    def __repr__(self):
        return "%s(%r, %s)" % (self.__class__.__name__, self.path, self.temporary)


class InputFile(_IOFile):
    pass


class TempInputFile(InputFile):
    def __init__(self, path: str):
        super().__init__(os.path.basename(path), temporary=True)


class OutputFile(_IOFile):
    pass


class TempOutputFile(OutputFile):
    def __init__(self, path: str):
        super().__init__(os.path.basename(path), temporary=True)


WrappedPipeType = Union[None, int, _IOFile, "AtomicCmd"]
PipeType = Union[str, WrappedPipeType]

OptionValueType = Union[str, float, "_IOFile", None]
OptionsType = Dict[
    str, Union[OptionValueType, List[OptionValueType], Tuple[OptionValueType, ...]]
]


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
        command: Iterable[Any],
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
            stdin -- A string path, InputFile, PIPE, DEVNULL, or an AtomicCmd instance.
                     If stdin is an AtomicCmd instance, stdin is read from its stdout.
            stdout -- A string path, InputFile, PIPE, DEVNULL, or an AtomicCmd instance.
            stderr -- A string path, InputFile, PIPE, DEVNULL, or an AtomicCmd instance.
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
        elif not isinstance(self._command[0], str):
            raise ValueError(self._command[0])
        elif self._command[0] not in self._executables:
            self._executables.add(self._command[0])

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
            if isinstance(value, _AtomicFile):
                self._record_atomic_file(value)
                if isinstance(value, Executable):
                    value = value.path
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

        errors = []  # type: List[str]
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

    def run(self, temp: str):
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

            # Explicitly set to DEVNULL to ensure that STDIN is not left open.
            if stdin is None:
                stdin = self.DEVNULL

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

    def join(self):
        """Similar to Popen.wait(), but returns the value wrapped in a list.
        Must be called before calling commit."""
        if not self._proc:
            return [None]

        self._running = False
        return_code = self._proc.wait()
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

    def commit(self, temp: str):
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

        committed_files = []  # type: List[str]
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
        if isinstance(value, Executable):
            return value.path
        elif self._set_cwd:
            if isinstance(value, InputFile):
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
            if isinstance(value, InputFile):
                if value.temporary:
                    return os.path.join(temp, value.path)

                return value.path
            elif isinstance(value, OutputFile):
                return fileutils.reroot_path(temp, value.path)
            elif isinstance(value, AuxilleryFile):
                return value.path
            elif value == "%(PYTHON)s":
                return sys.executable
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
        default: Optional[str] = None,
    ) -> WrappedPipeType:
        if pipe is None:
            if default is None:
                return None

            assert isinstance(self._command[0], str)
            executable = os.path.basename(self._command[0])
            # TODO: Use more sensible filename, e.g. log_{exec}_{counter}.{stdout/err}
            filename = "pipe_%s_%i.%s" % (executable, id(self), default)

            return filetype(filename, temporary=True)
        elif isinstance(pipe, int):
            return pipe
        elif isinstance(pipe, _AtomicFile):
            if isinstance(pipe, filetype):
                return pipe

            raise ValueError("expected %s, but got %s" % (filetype, pipe))
        elif isinstance(pipe, AtomicCmd):
            return pipe

        return filetype(pipe)

    @classmethod
    def _open_pipe(
        cls,
        temp_dir: str,
        pipe: WrappedPipeType,
        mode: str,
    ) -> Union[int, IO[bytes], None]:
        if pipe is None or isinstance(pipe, int):
            return pipe
        elif isinstance(pipe, AtomicCmd):
            return pipe._proc and pipe._proc.stdout

        if pipe.temporary or isinstance(pipe, OutputFile):
            return open(fileutils.reroot_path(temp_dir, pipe.path), mode)

        return open(pipe.path, mode)

    def __enter__(self):
        return self

    def __exit__(self, type: Any, _value: Any, _traceback: Any):
        self.terminate()
        self.join()

    def __str__(self):
        return atomicpp.pformat(self)


# The following ensures proper cleanup of child processes, for example in the
# case where multiprocessing.Pool.terminate() is called.
_PROCS = ()


def _cleanup_children(signum: int, _frame):
    for proc_ref in list(_PROCS):
        proc = proc_ref()
        if proc:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except OSError:
                # Ignore already closed processes, etc.
                pass
    sys.exit(-signum)


def _add_to_killlist(proc):
    global _PROCS

    if isinstance(_PROCS, tuple):
        signal.signal(signal.SIGTERM, _cleanup_children)
        _PROCS = set()

    _PROCS.add(weakref.ref(proc, _PROCS.remove))
