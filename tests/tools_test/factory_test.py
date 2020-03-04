#!/usr/bin/python
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
def test_paleomix_command():
    stdout, stderr = check_run(["paleomix"])

    assert "PALEOMIX - pipelines and tools for NGS data analyses" in stdout
    assert stderr == ""


FACTORY_COMMANDS = (
    (
        "bam_pipeline",
        "usage: paleomix bam_pipeline run [-h] [--version] [--log-file LOG_FILE]",
    ),
    (
        "trim_pipeline",
        "usage: paleomix trim_pipeline run [-h] [--version] [--log-file LOG_FILE]",
    ),
    (
        "phylo_pipeline",
        "usage: paleomix phylo_pipeline [-h] [--version] [--log-file LOG_FILE]",
    ),
    (
        "cleanup",
        "usage: paleomix cleanup --temp-prefix prefix --fasta reference.fasta < in.sam",
    ),
    ("coverage", "usage: paleomix coverage [options] sorted.bam [out.coverage]"),
    ("depths", "usage: paleomix depths [options] sorted.bam [out.depths]"),
    (
        "rmdup_collapsed",
        "usage: paleomix rmdup_collapsed [options] < sorted.bam > out.bam",
    ),
    (
        "gtf_to_bed",
        "usage: paleomix gtf_to_bed [options] in.gtf out_prefix [in.scaffolds]",
    ),
    ("vcf_filter", "usage: paleomix vcf_filter [-h] [--version] [--reset-filter]"),
    (
        "vcf_to_fasta",
        "usage: paleomix vcf_to_fasta [options] --genotype in.vcf --intervals in.bed",
    ),
)


# Simple test that all commands can be executed
@pytest.mark.parametrize("command, expected", FACTORY_COMMANDS)
def test_factory__commands(command, expected):
    cmd = factory.new(command)
    call = cmd.finalized_call
    if command in ("bam_pipeline", "trim_pipeline"):
        call.append("run")

    stdout, stderr = check_run(call + ["--help"])

    assert stdout.split("\n")[0] == expected
    assert stderr == ""
