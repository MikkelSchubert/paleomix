from atomiccmd import AtomicCmd, CmdError


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
