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
import os
import signal

import pytest

from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import ParallelCmds, SequentialCmds
from paleomix.atomiccmd.pprint import pformat, _pformat_list


###############################################################################
###############################################################################
# pformat


def test_pformat__simple():
    cmd = AtomicCmd(("touch", "something"))
    assert pformat(cmd) == (
        "Command = touch something\n"
        "STDOUT* = '${TEMP_DIR}/pipe_touch_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_touch_%i.stderr'"
    ) % (id(cmd), id(cmd))


def test_pformat__simple__running(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    assert pformat(cmd) == (
        "Command = sleep 10\n"
        "Status  = Running\n"
        "STDOUT* = '{temp_dir}/pipe_sleep_{id}.stdout'\n"
        "STDERR* = '{temp_dir}/pipe_sleep_{id}.stderr'\n"
        "CWD     = '{cwd}'"
    ).format(id=id(cmd), cwd=os.getcwd(), temp_dir=tmp_path)
    cmd.terminate()
    cmd.join()


def test_pformat__simple__running__set_cwd(tmp_path):
    cmd = AtomicCmd(("sleep", "10"), set_cwd=True)
    cmd.run(tmp_path)
    assert pformat(cmd) == (
        "Command = sleep 10\n"
        "Status  = Running\n"
        "STDOUT* = 'pipe_sleep_{id}.stdout'\n"
        "STDERR* = 'pipe_sleep_{id}.stderr'\n"
        "CWD     = '{temp_dir}'"
    ).format(id=id(cmd), temp_dir=tmp_path)
    cmd.terminate()
    cmd.join()


def test_pformat__simple__done(tmp_path):
    cmd = AtomicCmd("true")
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert pformat(cmd) == (
        "Command = true\n"
        "Status  = Exited with return-code 0\n"
        "STDOUT* = '{temp_dir}/pipe_true_{id}.stdout'\n"
        "STDERR* = '{temp_dir}/pipe_true_{id}.stderr'\n"
        "CWD     = '{cwd}'"
    ).format(id=id(cmd), cwd=os.getcwd(), temp_dir=tmp_path)


def test_pformat__simple__done__before_join(tmp_path):
    cmd = AtomicCmd("true")
    cmd.run(tmp_path)
    cmd._proc.wait()
    assert pformat(cmd) == (
        "Command = true\n"
        "Status  = Exited with return-code 0\n"
        "STDOUT* = '{temp_dir}/pipe_true_{id}.stdout'\n"
        "STDERR* = '{temp_dir}/pipe_true_{id}.stderr'\n"
        "CWD     = '{cwd}'"
    ).format(id=id(cmd), cwd=os.getcwd(), temp_dir=tmp_path)
    assert cmd.join() == [0]


def test_pformat__simple__done__set_cwd(tmp_path):
    cmd = AtomicCmd("true", set_cwd=True)
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert pformat(cmd) == (
        "Command = true\n"
        "Status  = Exited with return-code 0\n"
        "STDOUT* = 'pipe_true_{id}.stdout'\n"
        "STDERR* = 'pipe_true_{id}.stderr'\n"
        "CWD     = '{temp_dir}'"
    ).format(id=id(cmd), temp_dir=tmp_path)


def test_pformat__simple__terminated_by_pipeline(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    cmd.terminate()
    assert cmd.join() == ["SIGTERM"]
    assert pformat(cmd) == (
        "Command = sleep 10\n"
        "Status  = Automatically terminated by PALEOMIX\n"
        "STDOUT* = '{temp_dir}/pipe_sleep_{id}.stdout'\n"
        "STDERR* = '{temp_dir}/pipe_sleep_{id}.stderr'\n"
        "CWD     = '{cwd}'"
    ).format(id=id(cmd), temp_dir=tmp_path, cwd=os.getcwd())


def test_pformat__simple__killed_by_signal(tmp_path):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(tmp_path)
    os.killpg(cmd._proc.pid, signal.SIGTERM)
    assert cmd.join() == ["SIGTERM"]
    assert pformat(cmd) == (
        "Command = sleep 10\n"
        "Status  = Terminated with signal SIGTERM\n"
        "STDOUT* = '{temp_dir}/pipe_sleep_{id}.stdout'\n"
        "STDERR* = '{temp_dir}/pipe_sleep_{id}.stderr'\n"
        "CWD     = '{cwd}'"
    ).format(id=id(cmd), temp_dir=tmp_path, cwd=os.getcwd())


def test_pformat__simple__temp_root_in_arguments(tmp_path):
    cmd = AtomicCmd(("echo", "${TEMP_DIR}"))
    cmd.run(tmp_path)
    assert cmd.join() == [0]
    assert pformat(cmd) == (
        "Command = echo '{temp_dir}'\n"
        "Status  = Exited with return-code 0\n"
        "STDOUT* = '{temp_dir}/pipe_echo_{id}.stdout'\n"
        "STDERR* = '{temp_dir}/pipe_echo_{id}.stderr'\n"
        "CWD     = '{cwd}'"
    ).format(id=id(cmd), temp_dir=tmp_path, cwd=os.getcwd())


###############################################################################
###############################################################################
# INFILE


def test_pformat__atomiccmd__simple_with_infile():
    cmd = AtomicCmd(("cat", "%(IN_SOMETHING)s"), IN_SOMETHING="/etc/fstab")
    assert pformat(cmd) == (
        "Command = cat /etc/fstab\n"
        "STDOUT* = '${TEMP_DIR}/pipe_cat_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_cat_%i.stderr'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_infile__set_cwd():
    cmd = AtomicCmd(
        ("cat", "%(IN_SOMETHING)s"), IN_SOMETHING="/etc/fstab", set_cwd=True
    )
    assert pformat(cmd) == (
        "Command = cat /etc/fstab\n"
        "STDOUT* = 'pipe_cat_%i.stdout'\n"
        "STDERR* = 'pipe_cat_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_temp_infile():
    cmd = AtomicCmd(("cat", "%(TEMP_IN_FILE)s"), TEMP_IN_FILE="infile.txt")

    assert pformat(cmd) == (
        "Command = cat '${TEMP_DIR}/infile.txt'\n"
        "STDOUT* = '${TEMP_DIR}/pipe_cat_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_cat_%i.stderr'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_temp_infile__set_cwd():
    cmd = AtomicCmd(
        ("zcat", "%(TEMP_IN_FILE)s"), TEMP_IN_FILE="infile.gz", set_cwd=True
    )

    assert pformat(cmd) == (
        "Command = zcat infile.gz\n"
        "STDOUT* = 'pipe_zcat_%i.stdout'\n"
        "STDERR* = 'pipe_zcat_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd), id(cmd))


###############################################################################
###############################################################################
# OUTFILE


def test_pformat__atomiccmd__simple_with_outfile():
    cmd = AtomicCmd(("touch", "%(OUT_RC)s"), OUT_RC="/etc/bashrc")
    assert pformat(cmd) == (
        "Command = touch '${TEMP_DIR}/bashrc'\n"
        "STDOUT* = '${TEMP_DIR}/pipe_touch_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_touch_%i.stderr'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_outfile__set_cwd():
    cmd = AtomicCmd(("touch", "%(OUT_RC)s"), OUT_RC="/etc/bashrc", set_cwd=True)

    assert pformat(cmd) == (
        "Command = touch bashrc\n"
        "STDOUT* = 'pipe_touch_%i.stdout'\n"
        "STDERR* = 'pipe_touch_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_temp_outfile():
    cmd = AtomicCmd(("touch", "%(TEMP_OUT_RC)s"), TEMP_OUT_RC="bashrc")

    assert pformat(cmd) == (
        "Command = touch '${TEMP_DIR}/bashrc'\n"
        "STDOUT* = '${TEMP_DIR}/pipe_touch_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_touch_%i.stderr'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_temp_outfile__set_cwd():
    cmd = AtomicCmd(("touch", "%(TEMP_OUT_RC)s"), TEMP_OUT_RC="bashrc", set_cwd=True)

    assert pformat(cmd) == (
        "Command = touch bashrc\n"
        "STDOUT* = 'pipe_touch_%i.stdout'\n"
        "STDERR* = 'pipe_touch_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd), id(cmd))


###############################################################################
###############################################################################
# STDIN


def test_pformat__atomiccmd__simple_with_stdin():
    cmd = AtomicCmd("gzip", IN_STDIN="/etc/fstab")
    assert pformat(cmd) == (
        "Command = gzip\n"
        "STDIN   = '/etc/fstab'\n"
        "STDOUT* = '${TEMP_DIR}/pipe_gzip_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_gzip_%i.stderr'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_stdin__set_cwd():
    cmd = AtomicCmd("gzip", IN_STDIN="/etc/fstab", set_cwd=True)
    assert pformat(cmd) == (
        "Command = gzip\n"
        "STDIN   = '/etc/fstab'\n"
        "STDOUT* = 'pipe_gzip_%i.stdout'\n"
        "STDERR* = 'pipe_gzip_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_temp_stdin():
    cmd = AtomicCmd("gzip", TEMP_IN_STDIN="stabstabstab")
    assert pformat(cmd) == (
        "Command = gzip\n"
        "STDIN*  = '${TEMP_DIR}/stabstabstab'\n"
        "STDOUT* = '${TEMP_DIR}/pipe_gzip_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_gzip_%i.stderr'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_temp_stdin__set_cwd():
    cmd = AtomicCmd("gzip", TEMP_IN_STDIN="stabstabstab", set_cwd=True)
    assert pformat(cmd) == (
        "Command = gzip\n"
        "STDIN*  = 'stabstabstab'\n"
        "STDOUT* = 'pipe_gzip_%i.stdout'\n"
        "STDERR* = 'pipe_gzip_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd), id(cmd))


def test_pformat__atomiccmd__simple_with_stdin__cmd():
    cmd_1 = AtomicCmd("gzip", OUT_STDOUT=AtomicCmd.PIPE)
    cmd_2 = AtomicCmd("gzip", IN_STDIN=cmd_1)
    assert pformat(cmd_2) == (
        "Command = gzip\n"
        "STDIN   = <PIPE>\n"
        "STDOUT* = '${TEMP_DIR}/pipe_gzip_%i.stdout'\n"
        "STDERR* = '${TEMP_DIR}/pipe_gzip_%i.stderr'"
    ) % (id(cmd_2), id(cmd_2))


###############################################################################
###############################################################################
# STDOUT


def test_pformat__atomiccmd__simple_with_stdout():
    cmd = AtomicCmd(("echo", "Water. Water."), OUT_STDOUT="/dev/ls")
    assert pformat(cmd) == (
        "Command = echo 'Water. Water.'\n"
        "STDOUT  = '${TEMP_DIR}/ls'\n"
        "STDERR* = '${TEMP_DIR}/pipe_echo_%i.stderr'"
    ) % (id(cmd),)


def test_pformat__atomiccmd__simple_with_stdout__set_cwd():
    cmd = AtomicCmd(("echo", "*pant*. *pant*."), OUT_STDOUT="/dev/barf", set_cwd=True)
    assert pformat(cmd) == (
        "Command = echo '*pant*. *pant*.'\n"
        "STDOUT  = 'barf'\n"
        "STDERR* = 'pipe_echo_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd),)


def test_pformat__atomiccmd__simple_with_temp_stdout():
    cmd = AtomicCmd(("echo", "Oil. Oil."), TEMP_OUT_STDOUT="dm")
    assert pformat(cmd) == (
        "Command = echo 'Oil. Oil.'\n"
        "STDOUT* = '${TEMP_DIR}/dm'\n"
        "STDERR* = '${TEMP_DIR}/pipe_echo_%i.stderr'"
    ) % (id(cmd),)


def test_pformat__atomiccmd__simple_with_temp_stdout__set_cwd():
    cmd = AtomicCmd(
        ("echo", "Room service. Room service."), TEMP_OUT_STDOUT="pv", set_cwd=True
    )
    assert pformat(cmd) == (
        "Command = echo 'Room service. Room service.'\n"
        "STDOUT* = 'pv'\n"
        "STDERR* = 'pipe_echo_%i.stderr'\n"
        "CWD     = '${TEMP_DIR}'"
    ) % (id(cmd),)


def test_pformat__atomiccmd__simple_with_stdout_pipe():
    cmd = AtomicCmd(("echo", "!"), OUT_STDOUT=AtomicCmd.PIPE)
    assert pformat(cmd) == (
        "Command = echo '!'\n"
        "STDOUT  = <PIPE>\n"
        "STDERR* = '${TEMP_DIR}/pipe_echo_%i.stderr'"
    ) % (id(cmd),)


###############################################################################
###############################################################################
# ParallelCmds


@pytest.mark.parametrize(
    "cls, description",
    ((ParallelCmds, "Parallel processes"), (SequentialCmds, "Sequential processes")),
)
def test_pformat__sets__simple(cls, description):
    template = (
        "{description}:\n"
        "  Process 1:\n"
        "    Command = echo foo\n"
        "    STDOUT  = Piped to process 2\n"
        "    STDERR* = '${{TEMP_DIR}}/pipe_echo_{cmd_1_id}.stderr'\n"
        "\n"
        "  Process 2:\n"
        "    Command = gzip\n"
        "    STDIN   = Piped from process 1\n"
        "    STDOUT* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2_id}.stdout'\n"
        "    STDERR* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2_id}.stderr'"
    )

    cmd_1 = AtomicCmd(("echo", "foo"), OUT_STDOUT=AtomicCmd.PIPE)
    cmd_2 = AtomicCmd("gzip", IN_STDIN=cmd_1)
    cmd = cls((cmd_1, cmd_2))
    assert pformat(cmd) == template.format(
        description=description, cmd_1_id=id(cmd_1), cmd_2_id=id(cmd_2)
    )


def test_pformat__sets__nested():
    cmd_1 = AtomicCmd(("echo", "foo"), OUT_STDOUT=AtomicCmd.PIPE)
    cmd_2 = AtomicCmd("gzip", IN_STDIN=cmd_1)
    cmd_3 = AtomicCmd("sha1sum")
    set_1 = ParallelCmds((cmd_1, cmd_2))
    set_2 = SequentialCmds((set_1, cmd_3))
    assert pformat(set_2) == (
        "Sequential processes:\n"
        "  Parallel processes:\n"
        "    Process 1:\n"
        "      Command = echo foo\n"
        "      STDOUT  = Piped to process 2\n"
        "      STDERR* = '${{TEMP_DIR}}/pipe_echo_{cmd_1}.stderr'\n"
        "\n"
        "    Process 2:\n"
        "      Command = gzip\n"
        "      STDIN   = Piped from process 1\n"
        "      STDOUT* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2}.stdout'\n"
        "      STDERR* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2}.stderr'\n"
        "\n"
        "  Process 3:\n"
        "    Command = sha1sum\n"
        "    STDOUT* = '${{TEMP_DIR}}/pipe_sha1sum_{cmd_3}.stdout'\n"
        "    STDERR* = '${{TEMP_DIR}}/pipe_sha1sum_{cmd_3}.stderr'"
    ).format(cmd_1=id(cmd_1), cmd_2=id(cmd_2), cmd_3=id(cmd_3))


###############################################################################
###############################################################################
# Bad input


@pytest.mark.parametrize("value", (1, {}, ""))
def test_pformat__bad_input(value):
    with pytest.raises(TypeError):
        pformat(value)


###############################################################################
###############################################################################
# _pformat_list


def test_pformat_list__empty():
    assert _pformat_list([]) == ""


def test_pformat_list__single():
    assert _pformat_list([3]) == "3"


def test_pformat_list__multiple():
    assert _pformat_list([3, 2, 1]) == "3 2 1"


def test_pformat_list__wrapped():
    assert _pformat_list([3, 2, 1], width=1) == "3 \\\n    2 \\\n    1"
    assert _pformat_list([3, 2, 1], width=2) == "3 \\\n    2 \\\n    1"
    assert _pformat_list([3, 2, 1], width=3) == "3 \\\n    2 \\\n    1"
    assert _pformat_list([3, 2, 1], width=4) == "3 2 \\\n    1"
    assert _pformat_list([3, 2, 1], width=5) == "3 2 \\\n    1"
    assert _pformat_list([3, 2, 1], width=6) == "3 2 1"
    assert _pformat_list([3, 2, 1], width=7) == "3 2 1"


def test_pformat_list__escaped():
    assert _pformat_list(["a", "b c"], width=100) == "a 'b c'"
    assert _pformat_list(["a", "$c"], width=100) == "a '$c'"
    assert _pformat_list(["!a", "c"], width=100) == "'!a' c"
    assert _pformat_list(["a", "'c"], width=100) == """a ''"'"'c'"""
