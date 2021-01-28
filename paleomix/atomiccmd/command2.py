#!/usr/bin/env python3
import os
import signal

import paleomix.atomiccmd.pprint as atomicpp
import paleomix.common.fileutils as fileutils
import paleomix.common.procs as procs

from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.atomiccmd.command import CmdError, _add_to_killlist


class _AtomicFile:
    def __init__(self, path):
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
    def __init__(self, path, temporary=False):
        super().__init__(path)
        self.temporary = bool(temporary)

        if temporary and os.path.dirname(self.path):
            raise ValueError(self.path)

    def __repr__(self):
        return "%s(%r, %s)" % (self.__class__.__name__, self.path, self.temporary)


class InputFile(_IOFile):
    pass


class OutputFile(_IOFile):
    pass


class AtomicCmd2:
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

    PIPE = procs.PIPE
    DEVNULL = procs.DEVNULL

    InFile = InputFile
    OutFile = OutputFile
    AuxFile = AuxilleryFile
    Executable = Executable

    def __init__(
        self,
        command,
        stdin=None,
        stdout=None,
        stderr=None,
        set_cwd=False,
        extra_files=(),
        requirements=(),
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
            requirements -- An optional sequence of callables that can be used to check
                            that requirements are met for invoking this command.

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
        self._requirements = frozenset(requirements)
        self._input_files = set()
        self._output_files = set()
        self._output_files_map = {}
        self._auxiliary_files = set()

        self.append(*safe_coerce_to_tuple(command))
        if not self._command or not self._command[0]:
            raise ValueError("Empty command in AtomicCmd2 constructor")
        elif self._command[0] not in self._executables:
            self._executables.add(self._command[0])

        self._stdin = stdin = self._wrap_pipe(InputFile, stdin)
        self._stdout = stdout = self._wrap_pipe(OutputFile, stdout, "stdout")
        self._stderr = stderr = self._wrap_pipe(OutputFile, stderr, "stderr")

        pipes = []
        for pipe in (stdin, stdout, stderr):
            if isinstance(pipe, _AtomicFile):
                pipes.append(pipe)

        self.add_extra_files(extra_files)
        self.add_extra_files(pipes)

        for value in self._requirements:
            if not callable(value):
                raise TypeError("requirement must be callable, not %r" % (value,))

    def append(self, *args):
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

    def add_extra_files(self, files):
        if self._proc is not None:
            raise CmdError("cannot modify already started command")

        for value in files:
            if not isinstance(value, _AtomicFile):
                raise ValueError(value)

            self._record_atomic_file(value)

    def append_options(self, options, pred=lambda s: s.startswith("-")):
        if isinstance(options, dict):
            options = options.items()

        for (key, values) in options:
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

    @property
    def output_files(self):
        return frozenset(
            afile.path for afile in self._output_files if not afile.temporary
        )

    @property
    def expected_output_files(self):
        return frozenset(
            afile.basename() for afile in self._output_files if not afile.temporary
        )

    @property
    def optional_output_files(self):
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

    def run(self, temp):
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
            temp = "" if self._set_cwd else os.path.abspath(temp)
            call = self.to_call(temp)

            # Explicitly set to DEVNULL to ensure that STDIN is not left open.
            if stdin is None:
                stdin = self.DEVNULL

            self._proc = procs.open_proc(
                call,
                stdin=stdin,
                stdout=stdout,
                stderr=stderr,
                cwd=cwd,
                preexec_fn=os.setsid,
            )
        except Exception as error:
            message = "Error running commands:\n  Call = %r\n  Error = %r"
            raise CmdError(message % (self._command, error))
        finally:
            # Close pipes to allow the command to recieve SIGPIPE
            for handle in (stdin, stdout, stderr):
                if handle not in (None, self.PIPE, self.DEVNULL):
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

    def commit(self, temp):
        temp = fileutils.fspath(temp)
        if not self.ready():
            raise CmdError("Attempting to commit before command has completed")
        elif self._running:
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

        committed_files = []
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

    def to_call(self, temp):
        return [self._to_path(temp, value) for value in self._command]

    def _to_path(self, temp, value):
        if isinstance(value, InputFile):
            if value.temporary:
                return os.path.join(temp, value.path)
            elif self._set_cwd:
                return os.path.abspath(value.path)
            else:
                return value.path
        elif isinstance(value, OutputFile):
            if self._set_cwd:
                return os.path.basename(value.path)
            else:
                return fileutils.reroot_path(temp, value.path)
        else:
            return value.replace("%(TEMP_DIR)s", temp)

    def _record_atomic_file(self, value):
        if isinstance(value, AuxilleryFile):
            self._auxiliary_files.add(value.path)
        elif isinstance(value, Executable):
            self._executables.add(value.path)
        elif isinstance(value, InputFile):
            self._input_files.add(value)
        elif isinstance(value, OutputFile):
            basename = value.basename()
            # All output is saved in the same folder, so it is not possible to have
            # different output files with the same basename. Different OutputFile
            # objects with the same path are assume to refer to the same file
            if self._output_files_map.get(basename, value.path) != value.path:
                raise CmdError("multiple output files with name %r" % (basename,))

            self._output_files.add(value)
            self._output_files_map[basename] = value.path
        else:
            raise ValueError(value)

    def _wrap_pipe(self, filetype, pipe, default=None):
        if pipe is None:
            if default is None:
                return None

            executable = self._command[0]
            # TODO: Use more sensible filename, e.g. log_{exec}_{counter}.{stdout/err}
            filename = "pipe_%s_%i.%s" % (executable, id(self), default)

            return filetype(filename, temporary=True)
        elif pipe in (procs.PIPE, procs.DEVNULL):
            return pipe
        elif isinstance(pipe, _AtomicFile):
            if isinstance(pipe, filetype):
                return pipe

            raise ValueError("expected %s, but got %s" % (filetype, pipe))
        elif isinstance(pipe, AtomicCmd2):
            return pipe

        return filetype(pipe)

    @classmethod
    def _open_pipe(cls, temp_dir, pipe, mode):
        if pipe in (None, cls.PIPE, cls.DEVNULL):
            return pipe
        elif isinstance(pipe, AtomicCmd2):
            return pipe._proc and pipe._proc.stdout

        return open(fileutils.reroot_path(temp_dir, pipe.path), mode)

    def __enter__(self):
        return self

    def __exit__(self, type, _value, _traceback):
        self.terminate()
        self.join()

    def __str__(self):
        return atomicpp.pformat(self)
