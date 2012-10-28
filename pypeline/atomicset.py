#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
#
import time

from pypeline.atomiccmd import AtomicCmd, CmdError
from pypeline.common.utilities import safe_coerce_to_tuple


class _CommandSet:
    def __init__(self, commands):
        self._commands = safe_coerce_to_tuple(commands)
        if not self._commands:
            raise CmdError("Empty list passed to command set")

    def ready(self):
        return all(cmd.ready() for cmd in self._commands)

    def commit(self, temp):
        for command in self._commands:
            command.commit(temp)

    @property
    def executables(self):
        files = []
        for command in self._commands:
            files.extend(command.executables)
        return files

    @property
    def input_files(self):
        files = []
        for command in self._commands:
            files.extend(command.input_files)
        return files

    @property
    def output_files(self):
        files = []
        for command in self._commands:
            files.extend(command.output_files)
        return files

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

    def join(self):
        commands = list(enumerate(self._commands))
        return_codes = [[None]] * len(commands)
        while commands:
            for (index, command) in commands:
                if command.ready():
                    return_codes[index] = command.join()
                    commands.remove((index, command))
                elif any(any(codes) for codes in return_codes):
                    command.terminate()
                    return_codes[index] = ["SIGTERM"]
                    commands.remove((index, command))
                
            time.sleep(1)

        return sum(return_codes, [])

    def terminate(self):
        for command in self._commands:
            command.terminate()




class SequentialCmds(_CommandSet):
    """This class wraps a set of AtomicCmds, running them sequentially.
    This class therefore corresponds a set of lines in a bash script, 
    each of which invokes a forground job. For example:
    $ bcftools view snps.bcf | bgzip > snps.vcf.bgz
    $ tabix snps.vcf.bgz

    The list of commands may include any type of command. Note that
    the run function only returns once each sub-command has completed.
    A command is only executed if the previous command in the sequence
    was succesfully completed."""

    def __init__(self, commands):
        _CommandSet.__init__(self, commands)

        for command in self._commands:
            if not isinstance(command, (AtomicCmd, _CommandSet)):
                raise CmdError("ParallelCmds must only contain AtomicCmds or other ParallelCmds!")
        self._ready = False


    def run(self, temp):
        for command in self._commands:
            command.run(temp)

            if any(command.join()):
                break

        self._ready = True

    
    def ready(self):
        return self._ready


    def join(self):
        return_codes = []
        for command in self._commands:
            return_codes.extend(command.join())
        return return_codes

