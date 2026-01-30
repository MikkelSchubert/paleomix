# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
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
            f"Command returned {proc.returncode}: {call!r}:\n"
            f"STDOUT: {stdout!r}\n"
            f"STDERR: {stderr!r}"
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
