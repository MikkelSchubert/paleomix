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
import os
import subprocess

from nose.tools import \
     assert_in, \
     assert_equal

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
        raise ProcError("Command returned %i: %r:\nSTDOUT: %r\nSTDERR: %r"
                        % (proc.returncode, call, stdout, stderr))

    return stdout, stderr


# Simple test of the paleomxi command
def test_paleomix_command():
    stdout, stderr = check_run(["paleomix"])

    assert_equal("", stdout)
    assert_in("PALEOMIX - pipelines and tools for NGS data analyses.", stderr)


# Simple test that all commands can be executed
def test_factory__commands():
    def _do_test_factory__commands(command, expected):
        cmd = factory.new(command)
        call = cmd.finalized_call
        if command in ("bam_pipeline", "trim_pipeline"):
            call.append("run")

        stdout, stderr = check_run(call + ["--help"])

        assert_equal(expected, stdout.split("\n")[0])
        assert_equal("", stderr)

    commands = (("bam_pipeline",    "Usage: paleomix bam_pipeline <command> [options] [makefiles]"),
                ("trim_pipeline",   "Usage: paleomix trim_pipeline <command> [options] [makefiles]"),
                ("phylo_pipeline",  "Usage: paleomix phylo_pipeline <command> [options] [makefiles]"),
                ("cleanup",         "usage: paleomix cleanup --temp-prefix prefix --fasta reference.fasta < in.sam"),
                ("coverage",        "usage: paleomix coverage [options] sorted.bam [out.coverage]"),
                ("depths",          "usage: paleomix depths [options] sorted.bam [out.depths]"),
                ("duphist",         "usage: paleomix duphist sorted.bam > out.histogram"),
                ("rmdup_collapsed", "usage: paleomix rmdup_collapsed [options] < sorted.bam > out.bam"),
                ("genotype",        "usage: paleomix genotype [options] sorted.bam out.vcf.bgz"),
                ("gtf_to_bed",      "usage: paleomix gtf_to_bed [options] in.gtf out_prefix [in.scaffolds]"),
                ("sample_pileup",   "usage: paleomix sample_pileup [options] --genotype in.vcf --intervals in.bed > out.fasta"),
                ("vcf_filter",      "Usage: paleomix vcf_filter [options] [in1.vcf, ...]"),
                ("vcf_to_fasta",    "usage: paleomix vcf_to_fasta [options] --genotype in.vcf --intervals in.bed"),
                ("cat",             "usage: paleomix cat [-h] [--output OUTPUT] file [file ...]"))

    for command, expected in commands:
        yield _do_test_factory__commands, command, expected


# Simple test of aliased commands
def test_factory__command_alias():
    def _do_test_factory__command_alias(alias, command):
        alias, command = [alias], [command]
        if alias == ["bam_pipeline"] or alias == ["trim_pipeline"]:
            alias.append("run")
            command.append("run")

        stdout_1, stderr_1 = check_run(alias + ["--help"])
        stdout_2, stderr_2 = check_run(["paleomix"] + command + ["--help"])

        assert_equal(stderr_1, stderr_1)
        assert_equal(stderr_2, stderr_2)

    commands = (("bam_pipeline", "bam_pipeline"),
                ("bam_rmdup_collapsed", "rmdup_collapsed"),
                ("conv_gtf_to_bed", "gtf_to_bed"),
                ("phylo_pipeline", "phylo_pipeline"),
                ("trim_pipeline", "trim_pipeline"))

    for alias, command in commands:
        yield _do_test_factory__command_alias, alias, command
