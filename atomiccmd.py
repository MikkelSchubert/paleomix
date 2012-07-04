import os
import subprocess

import fileutils


class CmdError(RuntimeError):
    def __init__(self, msg):
        RuntimeError.__init__(self, msg)

    

class AtomicCmd:
    """Executes a command, only moving resulting files to the destination
    directory if the command was succesful. This helps prevent the 
    accidential use of partial files in downstream analysis, and eases 
    restarting of a pipeline following errors (no cleanup)."""
    PIPE = subprocess.PIPE

    def __init__(self, command, stdin = None, stdout = None, stderr = None, **kwargs):
        self.proc  = None
        self._cmd  = [str(field) for field in command]
        self._temp = None
        
        # Fill the in_files and out_files dictionaries 
        self._in_files, self._out_files = {}, {}
        self._group_io_by_prefix(kwargs)

        self._stdin  = stdin
        self._stderr = stderr
        self._stdout = stdout
        self._handles = {}

        if type(stdin) is str:
            self._in_files["PIPE_STDIN"] = stdin
        if type(stdout) is str:
            self._out_files["PIPE_STDOUT"] = stdout
        if (type(stderr) is str) and (stdout != stderr):
            self._out_files["PIPE_STDERR"] = stderr

        self._temp_files = None


    def run(self, temp):
        """Runs the given command, saving files in the specified temp folder. To 
        move files to their final destination, call commit(). Note that in contexts
        where the *Cmds classes are used, this function may block."""
        stdin  = self._open_pipe(temp, self._stdin , "rb")
        stdout = self._open_pipe(temp, self._stdout, "wb")
        stderr = self._open_pipe(temp, self._stderr, "wb")

        self._temp = temp
        self._temp_files = self._generate_temp_filenames(root = temp)
        
        kwords = dict(self._temp_files)
        for (key, filename) in self._in_files.iteritems():
            if key.startswith("TEMP_IN_"):
                filename = os.path.join(temp, filename)
            kwords[key] = filename
        kwords["TEMP_DIR"] = temp

        cmd = [(field % kwords) for field in self._cmd]

        self.proc = subprocess.Popen(cmd, 
                                     stdin  = stdin,
                                     stdout = stdout,
                                     stderr = stderr)


    def wait(self):
        """Equivalent to Popen.wait(), but returns the value wrapped in a list."""
        return [self.proc.wait()]


    def poll(self):
        """Equivalent to Popen.poll(), but returns the value wrapped in a list."""
        return [self.proc.poll()]


    def executables(self):
        """Returns a list of executables required for the AtomicCmd."""
        return [self._cmd[0]]


    def input_files(self):
        """Returns a list of input files that are required by the AtomicCmd."""
        for (key, filename) in self._in_files.iteritems():
            if not key.startswith("TEMP_"):
                yield filename
        

    def output_files(self):
        """Checks that the expected output files have been generated."""
        for (key, filename) in self._out_files.iteritems():
            if not key.startswith("TEMP_"):
                yield filename


    def temp_files(self):
        return list(self._temp_files.itervalues())


    def commit(self):
        assert (self.poll() is not None)

        # Close any implictly opened pipes
        for (_, handle) in self._handles.itervalues():
            handle.close()

        for key in self._temp_files:
            temp_filename = self._temp_files[key]
            final_filename = self._out_files[key]

            if key.startswith("TEMP_OUT_"):
                os.remove(temp_filename)
                continue

            fileutils.move_file(temp_filename, final_filename)

        self.proc = None
        return True


    def __str__(self):
        temp = self._temp or "${TEMP}"
        kwords = self._generate_temp_filenames(root = temp)
        kwords.update(self._in_files)
        kwords["TEMP_DIR"] = temp
        
        def describe_pipe(pipe, prefix):
            if type(pipe) is str:
                return "%s '%s'" % (prefix, pipe)
            elif isinstance(pipe, AtomicCmd):
                return "%s [AtomicCmd]" % (prefix)
            elif pipe == AtomicCmd.PIPE:
                return "%s [PIPE]" % prefix
            else:
                return ""

        if self._stdout != self._stderr:
            out  = describe_pipe(self._stdout, " >")
            out += describe_pipe(self._stderr, " 2>")
        else:
            out  = describe_pipe(self._stdout, " &>")
        

        command = [(field % kwords) for field in self._cmd]
        return "<'%s'%s%s" % (" ".join(command), describe_pipe(self._stdin, " <"), out)


    def _group_io_by_prefix(self, io_kwords):
        for (key, value) in io_kwords.iteritems():
            if value is None:
                continue # Simplifies use of optional arguments
            elif (value is not None) and (type(value) not in (str, unicode)):
                raise RuntimeError("Invalid input file '%s' for '%s' is not a string: %s" \
                                     % (key, self.__class__.__name__, value))

            if key.startswith("IN_") or key.startswith("TEMP_IN_"):
                self._in_files[key] = value
            elif key.startswith("OUT_") or key.startswith("TEMP_OUT_"):
                self._out_files[key] = value
            else:
                if not (key.startswith("IN_") or key.startswith("OUT_")):
                    raise CmdError("Command contains unclassified argument: '%s' -> '%s'" \
                                   % (self.__class__.__name__, key))


    def _open_pipe(self, root, pipe, mode):
        if isinstance(pipe, AtomicCmd):
            assert mode == "rb"
            return pipe.proc.stdout
        elif not (type(pipe) is str):
            return pipe

        pipe = os.path.basename(pipe)
        if pipe not in self._handles:
            self._handles[pipe] = (mode, open(os.path.join(root, pipe), mode))

        pipe_mode, pipe = self._handles[pipe]
        if pipe_mode != mode:
            raise CmdError("Attempting to open pipe with different modes: '%s' -> '%s'" \
                               % (self, pipe))

        return pipe


    def _generate_temp_filenames(self, root):
        filenames = {}
        for (key, filename) in self._out_files.iteritems():
            filenames[key] = os.path.join(root, os.path.basename(filename))
        return filenames




class _CommandSet:
    def __init__(self, commands):
        self._commands = list(commands)

    def wait(self):
        return_codes = []
        for command in self._commands:
            return_codes.extend(command.wait())
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

    def temp_files(self):
        files = []
        for command in self._commands:
            files.extend(command.temp_files())
        return files

    def commit(self):
        result = True
        for command in self._commands:
            result &= command.commit()
        return result

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

            if any(command.wait()):
                break
