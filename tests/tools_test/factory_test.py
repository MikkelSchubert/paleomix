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

import subprocess

import pytest

import paleomix
import paleomix.__main__ as main
from paleomix.tools import factory

###############################################################################
###############################################################################


class ProcError(RuntimeError):
    pass


def check_run(call: list[str], expected_returncode: int = 0) -> tuple[str, str]:
    proc = subprocess.Popen(
        call,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.DEVNULL,
        encoding="utf-8",
        errors="replace",
        close_fds=True,
    )

    stdout, stderr = proc.communicate()
    if proc.returncode != expected_returncode:
        raise ProcError(
            "Command returned %i: %r:\nSTDOUT: %r\nSTDERR: %r"
            % (proc.returncode, call, stdout, stderr)
        )

    return stdout, stderr


# Simple test of the paleomxi command
@pytest.mark.slow
def test_paleomix_command() -> None:
    stdout, stderr = check_run(["paleomix"])

    assert "PALEOMIX - pipelines and tools for NGS data analyses" in stdout
    assert stderr == ""


# Simple test that all commands can be executed
@pytest.mark.slow
@pytest.mark.parametrize("command", main.COMMANDS)
def test_factory__command_usage(command: str) -> None:
    cmd = factory.new(command)
    call = cmd.to_call("%(TEMP_DIR)s")

    stdout, stderr = check_run([*call, "--help"])

    name = command.replace("_pipeline", "")
    assert stdout.startswith(f"usage: paleomix {name}")
    assert not stderr


# Simple test that all commands support '--version'
@pytest.mark.slow
@pytest.mark.parametrize("arg", ["-v", "--version"])
@pytest.mark.parametrize("command", main.COMMANDS)
def test_factory__command_versions(arg: str, command: str) -> None:
    cmd = factory.new(command)
    call = cmd.to_call("%(TEMP_DIR)s")

    stdout, stderr = check_run([*call, arg])

    name = command.replace("_pipeline", "")
    assert stdout.startswith(f"paleomix {name}")
    assert stdout.endswith(f" v{paleomix.__version__}\n")
    assert not stderr
