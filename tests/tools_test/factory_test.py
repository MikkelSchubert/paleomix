#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import os
import subprocess

import pytest

import paleomix
import paleomix.main as main
import paleomix.tools.factory as factory


###############################################################################
###############################################################################


class ProcError(RuntimeError):
    pass


def check_run(call, *args, **kwargs):
    devnull = os.open(os.devnull, os.O_RDONLY)
    kwargs.setdefault("stdin", devnull)
    kwargs.setdefault("close_fds", True)
    kwargs["stdout"] = subprocess.PIPE
    kwargs["stderr"] = subprocess.PIPE

    returncode = kwargs.pop("expected_returncode", 0)

    proc = subprocess.Popen(call, *args, **kwargs)
    os.close(devnull)

    stdout, stderr = proc.communicate()
    if proc.returncode != returncode:
        raise ProcError(
            "Command returned %i: %r:\nSTDOUT: %r\nSTDERR: %r"
            % (proc.returncode, call, stdout, stderr)
        )

    return stdout.decode("utf-8", "replace"), stderr.decode("utf-8", "replace")


# Simple test of the paleomxi command
@pytest.mark.slow
def test_paleomix_command():
    stdout, stderr = check_run(["paleomix"])

    assert "PALEOMIX - pipelines and tools for NGS data analyses" in stdout
    assert stderr == ""


# Simple test that all commands can be executed
@pytest.mark.slow
@pytest.mark.parametrize("command", main._COMMANDS)
def test_factory__command_usage(command):
    cmd = factory.new(command)
    call = cmd.finalized_call

    stdout, stderr = check_run(call + ["--help"])

    assert stdout.startswith("usage: paleomix {}".format(command))
    assert not stderr


# Simple test that all commands support '--version'
@pytest.mark.slow
@pytest.mark.parametrize("arg", ("-v", "--version"))
@pytest.mark.parametrize("command", main._COMMANDS)
def test_factory__command_versions(arg, command):
    cmd = factory.new(command)
    call = cmd.finalized_call

    stdout, stderr = check_run(call + [arg])

    assert stdout.startswith("paleomix {}".format(command))
    assert stdout.endswith(" v{}\n".format(paleomix.__version__))
    assert not stderr
