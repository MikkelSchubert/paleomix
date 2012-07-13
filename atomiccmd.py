import os
import sys
import signal
import types
import weakref
import subprocess

import fileutils



_PIPES = ('STDIN', 'STDOUT', 'STDERR')
_PREFIXES = ('IN_', 'TEMP_IN', 'OUT_', 'TEMP_OUT')




class CmdError(RuntimeError):
    def __init__(self, msg):
        RuntimeError.__init__(self, msg)

    

class AtomicCmd:
    """Executes a command, only moving resulting files to the destination
    directory if the command was succesful. This helps prevent the 
    accidential use of partial files in downstream analysis, and eases 
    restarting of a pipeline following errors (no cleanup).

    When an AtomicCmd is run(), a signal handler is installed for SIGTERM,
    which ensures that any running processes are terminated. In the absence
    of this, AtomicCmds run in terminated subprocesses can result in still
    running children after the termination of the parents."""
    PIPE = subprocess.PIPE

    def __init__(self, command, stdin = None, stdout = None, stderr = None, **kwargs):
        """
        
        If the above is prefixed with "TEMP_", the files are read from / written
        to the temporary folder in which the command is executed. Note that all
        TEMP_OUT_ files are deleted when commit is called (if they exist).

        In addition, the follow special names may be used with the above:
           STDIN  -- Takes a filename, or an AtomicCmd, in which case the stdout
                     of that command is piped to the stdin of this instance.
           STDOUT -- Takes a filename, or the special value PIPE to allow
                     another AtomicCmd instance to use the output directly.
           STDERR -- Takes a filename. If the filename is the same as STDOUT,
                     all output is redirected to that file.
                     
        Each pipe can only be used once (either OUT_ or TEMP_OUT_)."""

        self._proc    = None
        self._command = [str(field) for field in command]
        self._handles = {}
        self._files   = {}
        
        # FIXME: Drop stdin, stderr, stdout parameters from constructor
        kwargs["IN_STDIN"]  = stdin
        kwargs["OUT_STDOUT"] = stdout
        kwargs["OUT_STDERR"] = stderr
        self._process_arguments(kwargs)

        

    def run(self, temp):
        """Runs the given command, saving files in the specified temp folder. To 
        move files to their final destination, call commit(). Note that in contexts
        where the *Cmds classes are used, this function may block."""

        kwords = self._generate_filenames(root = temp)
        stdin  = self._open_pipe(kwords, "IN_STDIN" , "rb")
        stdout = self._open_pipe(kwords, "OUT_STDOUT", "wb")
        stderr = self._open_pipe(kwords, "OUT_STDERR", "wb")
        command = [(field % kwords) for field in self._command]

        self._proc = subprocess.Popen(command, 
                                      stdin  = stdin,
                                      stdout = stdout,
                                      stderr = stderr)
        
        # Allow subprocesses to be killed in case of a SIGTERM
        _add_to_killlist(self._proc)


    def join(self):
        """Similar to Popen.wait(), but returns the value wrapped in a list,
        and ensures that any opened handles are closed. Must be called before
        calling commit."""
        try:
            return_codes = [self._proc.wait()]
        finally:
            # Close any implictly opened pipes
            for (mode, handle) in self._handles.itervalues():
                if "w" in mode:
                    handle.flush()
                handle.close()
            self._handles = {}

        return return_codes


    def poll(self):
        """Equivalent to Popen.poll(), but returns the value wrapped in a list."""
        return [self._proc.poll()]


    def executables(self):
        """Returns a list of executables required for the AtomicCmd."""
        return [self._command[0]]


    def input_files(self):
        """Returns a list of input files that are required by the AtomicCmd."""
        for (key, filename) in self._files.iteritems():
            if key.startswith("IN_") and isinstance(filename, types.StringTypes):
                yield filename
        

    def output_files(self):
        """Checks that the expected output files have been generated."""
        for (key, filename) in self._files.iteritems():
            if key.startswith("OUT_") and isinstance(filename, types.StringTypes):
                yield filename


    def commit(self, temp):
        assert (self.poll() is not None)
        assert not self._handles

        self._proc = None

        temp_files = self._generate_filenames(root = temp)
        for (key, filename) in temp_files.iteritems():
            if isinstance(filename, types.StringTypes):
                if key.startswith("OUT_") and not os.path.exists(filename):
                    raise CmdError("Command did not create expected output file: " + filename)

        for (key, filename) in temp_files.iteritems():
            if isinstance(filename, types.StringTypes):
                if key.startswith("OUT_"):
                    fileutils.move_file(filename, self._files[key])
                elif key.startswith("TEMP_OUT_") and os.path.exists(filename):
                    os.remove(filename)


    def __str__(self):        
        def describe_pipe(template, pipe):
            if isinstance(pipe, types.StringTypes):
                return template % pipe
            elif isinstance(pipe, AtomicCmd):
                return template % "[AtomicCmd]"
            elif (pipe == AtomicCmd.PIPE):
                return template % "[PIPE]"
            else:
                return ""

        kwords = self._generate_filenames(root = "${TEMP}")
        command = " ".join([(field % kwords) for field in self._command])

        stdin = self._files.get("IN_STDIN")
        if stdin:
            command += describe_pipe(" < %s", stdin)

        stdout = self._files.get("OUT_STDOUT")
        stderr = self._files.get("OUT_STDERR")
        if (stdout != stderr):
            if stdout:
                command += describe_pipe(" > %s", stdout)
            if stderr:
                command += describe_pipe(" 2> %s", stderr)
        elif stdout:
            command += describe_pipe(" &> %s", stdout)
       
        return "<%s>" % command


    def _process_arguments(self, kwargs):
        self._validate_pipes(kwargs)

        for (key, value) in kwargs.iteritems():
            if self._validate_argument(key, value):
                self._files[key] = value

        executable = os.path.basename(self._command[0])
        for pipe in ("STDOUT", "STDERR"):
            if not (kwargs.get("OUT_" + pipe) or kwargs.get("TEMP_OUT_" + pipe)):
                filename = "pipe_%s_%i.%s" % (executable, id(self), pipe.lower())
                self._files["TEMP_OUT_" + pipe] = filename


    @classmethod
    def _validate_pipes(cls, kwargs):
        """Checks that no single pipe is specified multiple times, e.i. being specified
        both for a temporary and a final (outside the temp dir) file. For example,
        either IN_STDIN or TEMP_IN_STDIN must be specified, but not both."""
        for pipe in _PIPES:
            if (kwargs.get("IN_" + pipe) and kwargs.get("TEMP_IN_" + pipe)):
                raise CmdError, "Pipe (%s) must be specified at most once." % pipe


    @classmethod
    def _validate_argument(cls, key, value):
        if (value == cls.PIPE) and key not in ("OUT_STDOUT", "TEMP_OUT_STDOUT"):
            raise ValueError, "PIPE is only allow for *_STDOUT, not " + key
        elif not (value is None or isinstance(value, types.StringTypes)):
            if not ((key == "IN_STDIN") and isinstance(value, AtomicCmd)):
                ValueError("Invalid input file '%s' for '%s' is not a string: %s" \
                                 % (key, cls.__name__, value))
        elif not any(key.startswith(prefix) for prefix in _PREFIXES):
            raise CmdError("Command contains invalid argument: '%s' -> '%s'" \
                               % (cls.__name__, key))
        
        # Values are allowed to be none, but such entries are skipped
        # This simplifies the use of keyword parameters with no default values.
        return bool(value)


    def _open_pipe(self, kwords, pipe, mode):
        filename = kwords.get(pipe, kwords.get("TEMP_" + pipe))
        if filename in (None, self.PIPE):
            return filename
        elif isinstance(filename, AtomicCmd):
            return filename._proc.stdout
        
        if filename not in self._handles:
            self._handles[pipe] = (mode, open(filename, mode))

        pipe_mode, pipe = self._handles[pipe]
        if pipe_mode != mode:
            raise CmdError("Attempting to open pipe with multiple modes: '%s' -> '%s'" \
                               % (self, filename))

        return pipe


    def _generate_filenames(self, root):
        filenames = {"TEMP_DIR" : root}
        for (key, filename) in self._files.iteritems():
            if isinstance(filename, types.StringTypes):
                if key.startswith("TEMP_") or key.startswith("OUT_"):
                    filename = os.path.join(root, os.path.basename(filename))
            filenames[key] = filename

        return filenames




class _CommandSet:
    def __init__(self, commands):
        self._commands = list(commands)

    def join(self):
        return_codes = []
        for command in self._commands:
            return_codes.extend(command.join())
        return return_codes

    def poll(self):
        return_codes = []
        for command in self._commands:
            return_codes.extend(command.poll())
        return return_codes

    def executables(self):
        files = []
        for command in self._commands:
            files.extend(command.executables())
        return files

    def input_files(self):
        files = []
        for command in self._commands:
            files.extend(command.input_files())
        return files

    def output_files(self):
        files = []
        for command in self._commands:
            files.extend(command.output_files())
        return files

    def commit(self, temp):
        for command in self._commands:
            command.commit(temp)

    def __str__(self):
        return "[%s]" % ", ".join(str(command) for command in self._commands)




class ParallelCmds(_CommandSet):
    """This class wraps a set of AtomicCmds, running them in parallel.
    This corresponds to a set of piped commands, which only terminate
    when all parts of the pipe have terminated. For example:
    $ dmesg | grep -i segfault | gzip > log.txt.gz

    Note that only AtomicCmds and ParallelCmds are allowed as 
    sub-commands for this class, since the model requires non-
    blocking commands."""

    def __init__(self, commands):
        _CommandSet.__init__(self, commands)

        for command in self._commands:
            if not isinstance(command, (AtomicCmd, ParallelCmds)):
                raise CmdError("ParallelCmds must only contain AtomicCmds or other ParallelCmds!")

    def run(self, temp):
        for command in self._commands:
            command.run(temp)




class SequentialCmds(_CommandSet):
    """This class wraps a set of AtomicCmds, running them sequentially.
    This class therefore corresponds a set of lines in a bash script, 
    each of which invokes a foreground task. For example:
    $ bcftools view snps.bcf | bgzip > snps.vcf.bgz
    $ tabix snps.vcf.bgz

    The list of commands may include any type of command. Note that
    the run function only returns once each sub-command has completed."""

    def __init__(self, commands):
        _CommandSet.__init__(self, commands)

        for command in self._commands:
            if not isinstance(command, (AtomicCmd, ParallelCmds)):
                raise CmdError("ParallelCmds must only contain AtomicCmds or other ParallelCmds!")


    def run(self, temp):
        for command in self._commands:
            command.run(temp)

            if any(command.join()):
                break



## The following ensures proper cleanup of child processes, for example in the case
## where multiprocessing.Pool.terminate() is called. The current implementation will
## leak weakref objects, but given the low number of commands called over the lifetime
## of a pipeline, this is considered acceptable for now. FIXME
_PROCS = []

def _cleanup_children(signum, _frame):
    for proc_ref in _PROCS:
        proc = proc_ref()
        if proc:
            proc.terminate()
    sys.exit(-signum)

def _add_to_killlist(proc):
    if not _PROCS:
        signal.signal(signal.SIGTERM, _cleanup_children)

    _PROCS.append(weakref.ref(proc))
