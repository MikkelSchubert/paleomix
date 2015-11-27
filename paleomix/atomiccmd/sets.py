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
import collections

import paleomix.atomiccmd.pprint as atomicpp

from paleomix.atomiccmd.command import AtomicCmd, CmdError
from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.common.fileutils import try_remove


class _CommandSet:
    def __init__(self, commands):
        self._commands = safe_coerce_to_tuple(commands)
        if not self._commands:
            raise CmdError("Empty list passed to command set")

        self._validate_commands()

    def commit(self, temp):
        committed_files = set()
        try:
            for command in self._commands:
                command.commit(temp)
                committed_files.update(command.output_files)
        except:
            # Cleanup after failed commit
            for fpath in committed_files:
                try_remove(fpath)
            raise

    def _collect_properties(key):  # pylint: disable=E0213
        def _collector(self):
            values = set()
            for command in self._commands:  # pylint: disable=W0212
                values.update(getattr(command, key))
            return values
        return property(_collector)

    input_files     = _collect_properties("input_files")
    output_files    = _collect_properties("output_files")
    auxiliary_files = _collect_properties("auxiliary_files")
    executables     = _collect_properties("executables")
    requirements    = _collect_properties("requirements")
    expected_temp_files = _collect_properties("expected_temp_files")
    optional_temp_files = _collect_properties("optional_temp_files")

    @property
    def stdout(self):
        raise CmdError("%s does not implement property 'stdout'!" \
                       % (self.__class__.__name__,))

    def terminate(self):
        for command in self._commands:
            command.terminate()

    def __str__(self):
        return atomicpp.pformat(self)

    def _validate_commands(self):
        if len(self._commands) != len(set(self._commands)):
            raise ValueError("Same command included multiple times in %s" \
                             % (self.__class__.__name__,))

        filenames = collections.defaultdict(int)
        for command in self._commands:
            for filename in command.expected_temp_files:
                filenames[filename] += 1
            for filename in command.optional_temp_files:
                filenames[filename] += 1

        clobbered = [filename for (filename, count) in filenames.items() if (count > 1)]
        if any(clobbered):
            raise CmdError("Commands clobber each others' files: %s" % (", ".join(clobbered),))


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

    def __init__(self, commands):
        self._joinable = False

        commands = safe_coerce_to_tuple(commands)
        for command in commands:
            if not isinstance(command, (AtomicCmd, ParallelCmds)):
                raise CmdError("ParallelCmds must only contain AtomicCmds or other ParallelCmds!")
        _CommandSet.__init__(self, commands)

    def run(self, temp):
        for command in self._commands:
            command.run(temp)
        self._joinable = True

    def ready(self):
        return all(cmd.ready() for cmd in self._commands)

    def join(self):
        sleep_time = 0.05
        commands = list(enumerate(self._commands))
        return_codes = [[None]] * len(commands)
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
        return sum(return_codes, [])




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

    def __init__(self, commands):
        self._ready = False

        commands = safe_coerce_to_tuple(commands)
        for command in commands:
            if not isinstance(command, (AtomicCmd, _CommandSet)):
                raise CmdError("ParallelCmds must only contain AtomicCmds or other ParallelCmds!")
        _CommandSet.__init__(self, commands)

    def run(self, temp):
        self._ready = False
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
