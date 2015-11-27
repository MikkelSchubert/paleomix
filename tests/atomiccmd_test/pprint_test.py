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
import StringIO

import nose
from nose.tools import assert_equal
from paleomix.common.testing import with_temp_folder

from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import ParallelCmds, SequentialCmds
from paleomix.atomiccmd.pprint import pprint, pformat, _pformat_list


###############################################################################
###############################################################################
# pprint

def test_pprint__simple():
    obj = StringIO.StringIO()
    cmd = AtomicCmd(("touch", "something"))
    pprint(cmd, out=obj)
    assert_equal(obj.getvalue(), ("<Command = ['touch', 'something']\n"
                                  " STDOUT* = '${TEMP_DIR}/pipe_touch_%i.stdout'\n"
                                  " STDERR* = '${TEMP_DIR}/pipe_touch_%i.stderr'>\n") % (id(cmd), id(cmd)))


###############################################################################
###############################################################################
# INFILE

def test_pformat__simple():
    cmd = AtomicCmd(("touch", "something"))
    assert_equal(pformat(cmd), ("<Command = ['touch', 'something']\n"
                                " STDOUT* = '${TEMP_DIR}/pipe_touch_%i.stdout'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_touch_%i.stderr'>") % (id(cmd), id(cmd)))


@with_temp_folder
def test_pformat__simple__running(temp_folder):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(temp_folder)
    assert_equal(pformat(cmd), ("<Command = ['sleep', '10']\n"
                                " Status  = Running ...\n"
                                " STDOUT* = '{temp_dir}/pipe_sleep_{id}.stdout'\n"
                                " STDERR* = '{temp_dir}/pipe_sleep_{id}.stderr'\n"
                                " CWD     = '{cwd}'>").format(id=id(cmd),
                                                              cwd=os.getcwd(),
                                                              temp_dir=temp_folder))
    cmd.terminate()
    cmd.join()


@with_temp_folder
def test_pformat__simple__running__set_cwd(temp_folder):
    cmd = AtomicCmd(("sleep", "10"), set_cwd=True)
    cmd.run(temp_folder)
    assert_equal(pformat(cmd), ("<Command = ['sleep', '10']\n"
                                " Status  = Running ...\n"
                                " STDOUT* = 'pipe_sleep_{id}.stdout'\n"
                                " STDERR* = 'pipe_sleep_{id}.stderr'\n"
                                " CWD     = '{temp_dir}'>").format(id=id(cmd),
                                                                   temp_dir=temp_folder))
    cmd.terminate()
    cmd.join()


@with_temp_folder
def test_pformat__simple__done(temp_folder):
    cmd = AtomicCmd("true")
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    assert_equal(pformat(cmd), ("<Command = ['true']\n"
                                " Status  = Exited with return-code 0\n"
                                " STDOUT* = '{temp_dir}/pipe_true_{id}.stdout'\n"
                                " STDERR* = '{temp_dir}/pipe_true_{id}.stderr'\n"
                                " CWD     = '{cwd}'>").format(id=id(cmd),
                                                              cwd=os.getcwd(),
                                                              temp_dir=temp_folder))


@with_temp_folder
def test_pformat__simple__done__before_join(temp_folder):
    cmd = AtomicCmd("true")
    cmd.run(temp_folder)
    cmd._proc.wait()
    assert_equal(pformat(cmd), ("<Command = ['true']\n"
                                " Status  = Exited with return-code 0\n"
                                " STDOUT* = '{temp_dir}/pipe_true_{id}.stdout'\n"
                                " STDERR* = '{temp_dir}/pipe_true_{id}.stderr'\n"
                                " CWD     = '{cwd}'>").format(id=id(cmd),
                                                              cwd=os.getcwd(),
                                                              temp_dir=temp_folder))
    assert_equal(cmd.join(), [0])


@with_temp_folder
def test_pformat__simple__done__set_cwd(temp_folder):
    cmd = AtomicCmd("true", set_cwd=True)
    cmd.run(temp_folder)
    assert_equal(cmd.join(), [0])
    assert_equal(pformat(cmd), ("<Command = ['true']\n"
                                " Status  = Exited with return-code 0\n"
                                " STDOUT* = 'pipe_true_{id}.stdout'\n"
                                " STDERR* = 'pipe_true_{id}.stderr'\n"
                                " CWD     = '{temp_dir}'>").format(id=id(cmd),
                                                                   temp_dir=temp_folder))


@with_temp_folder
def test_pformat__simple__killed(temp_folder):
    cmd = AtomicCmd(("sleep", "10"))
    cmd.run(temp_folder)
    cmd.terminate()
    assert_equal(cmd.join(), ["SIGTERM"])
    assert_equal(pformat(cmd), ("<Command = ['sleep', '10']\n"
                                " Status  = Terminated with signal SIGTERM\n"
                                " STDOUT* = '{temp_dir}/pipe_sleep_{id}.stdout'\n"
                                " STDERR* = '{temp_dir}/pipe_sleep_{id}.stderr'\n"
                                " CWD     = '{cwd}'>").format(id=id(cmd),
                                                              temp_dir=temp_folder,
                                                              cwd=os.getcwd()))


###############################################################################
###############################################################################
# INFILE

def test_pformat__atomiccmd__simple_with_infile():
    cmd = AtomicCmd(("cat", "%(IN_SOMETHING)s"), IN_SOMETHING="/etc/fstab")
    assert_equal(pformat(cmd), ("<Command = ['cat', '/etc/fstab']\n"
                                " STDOUT* = '${TEMP_DIR}/pipe_cat_%i.stdout'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_cat_%i.stderr'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_infile__set_cwd():
    cmd = AtomicCmd(("cat", "%(IN_SOMETHING)s"),
                    IN_SOMETHING="/etc/fstab",
                    set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['cat', '/etc/fstab']\n"
                                " STDOUT* = 'pipe_cat_%i.stdout'\n"
                                " STDERR* = 'pipe_cat_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_temp_infile():
    cmd = AtomicCmd(("cat", "%(TEMP_IN_FILE)s"),
                    TEMP_IN_FILE="infile.txt")
    assert_equal(pformat(cmd), ("<Command = ['cat', '${TEMP_DIR}/infile.txt']\n"
                                " STDOUT* = '${TEMP_DIR}/pipe_cat_%i.stdout'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_cat_%i.stderr'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_temp_infile__set_cwd():
    cmd = AtomicCmd(("zcat", "%(TEMP_IN_FILE)s"),
                    TEMP_IN_FILE="infile.gz",
                    set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['zcat', 'infile.gz']\n"
                                " STDOUT* = 'pipe_zcat_%i.stdout'\n"
                                " STDERR* = 'pipe_zcat_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd), id(cmd)))


###############################################################################
###############################################################################
# OUTFILE

def test_pformat__atomiccmd__simple_with_outfile():
    cmd = AtomicCmd(("touch", "%(OUT_RC)s"), OUT_RC="/etc/bashrc")
    assert_equal(pformat(cmd), ("<Command = ['touch', '${TEMP_DIR}/bashrc']\n"
                                " STDOUT* = '${TEMP_DIR}/pipe_touch_%i.stdout'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_touch_%i.stderr'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_outfile__set_cwd():
    cmd = AtomicCmd(("touch", "%(OUT_RC)s"), OUT_RC="/etc/bashrc", set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['touch', 'bashrc']\n"
                                " STDOUT* = 'pipe_touch_%i.stdout'\n"
                                " STDERR* = 'pipe_touch_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_temp_outfile():
    cmd = AtomicCmd(("touch", "%(TEMP_OUT_RC)s"), TEMP_OUT_RC="bashrc")
    assert_equal(pformat(cmd), ("<Command = ['touch', '${TEMP_DIR}/bashrc']\n"
                                " STDOUT* = '${TEMP_DIR}/pipe_touch_%i.stdout'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_touch_%i.stderr'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_temp_outfile__set_cwd():
    cmd = AtomicCmd(("touch", "%(TEMP_OUT_RC)s"), TEMP_OUT_RC="bashrc", set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['touch', 'bashrc']\n"
                                " STDOUT* = 'pipe_touch_%i.stdout'\n"
                                " STDERR* = 'pipe_touch_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd), id(cmd)))


###############################################################################
###############################################################################
# STDIN

def test_pformat__atomiccmd__simple_with_stdin():
    cmd = AtomicCmd("gzip", IN_STDIN="/etc/fstab")
    assert_equal(pformat(cmd), ("<Command = ['gzip']\n"
                                " STDIN   = '/etc/fstab'\n"
                                " STDOUT* = '${TEMP_DIR}/pipe_gzip_%i.stdout'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_gzip_%i.stderr'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_stdin__set_cwd():
    cmd = AtomicCmd("gzip", IN_STDIN="/etc/fstab", set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['gzip']\n"
                                " STDIN   = '/etc/fstab'\n"
                                " STDOUT* = 'pipe_gzip_%i.stdout'\n"
                                " STDERR* = 'pipe_gzip_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_temp_stdin():
    cmd = AtomicCmd("gzip", TEMP_IN_STDIN="stabstabstab")
    assert_equal(pformat(cmd), ("<Command = ['gzip']\n"
                                " STDIN*  = '${TEMP_DIR}/stabstabstab'\n"
                                " STDOUT* = '${TEMP_DIR}/pipe_gzip_%i.stdout'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_gzip_%i.stderr'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_temp_stdin__set_cwd():
    cmd = AtomicCmd("gzip", TEMP_IN_STDIN="stabstabstab", set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['gzip']\n"
                                " STDIN*  = 'stabstabstab'\n"
                                " STDOUT* = 'pipe_gzip_%i.stdout'\n"
                                " STDERR* = 'pipe_gzip_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd), id(cmd)))


def test_pformat__atomiccmd__simple_with_stdin__cmd():
    cmd_1 = AtomicCmd("gzip", OUT_STDOUT=AtomicCmd.PIPE)
    cmd_2 = AtomicCmd("gzip", IN_STDIN=cmd_1)
    assert_equal(pformat(cmd_2), ("<Command = ['gzip']\n"
                                  " STDIN   = <PIPE>\n"
                                  " STDOUT* = '${TEMP_DIR}/pipe_gzip_%i.stdout'\n"
                                  " STDERR* = '${TEMP_DIR}/pipe_gzip_%i.stderr'>") % (id(cmd_2), id(cmd_2)))


###############################################################################
###############################################################################
# STDOUT

def test_pformat__atomiccmd__simple_with_stdout():
    cmd = AtomicCmd(("echo", "Water. Water."), OUT_STDOUT="/dev/ls")
    assert_equal(pformat(cmd), ("<Command = ['echo', 'Water. Water.']\n"
                                " STDOUT  = '${TEMP_DIR}/ls'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_echo_%i.stderr'>") % (id(cmd),))


def test_pformat__atomiccmd__simple_with_stdout__set_cwd():
    cmd = AtomicCmd(("echo", "*pant*. *pant*."), OUT_STDOUT="/dev/barf", set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['echo', '*pant*. *pant*.']\n"
                                " STDOUT  = 'barf'\n"
                                " STDERR* = 'pipe_echo_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd),))


def test_pformat__atomiccmd__simple_with_temp_stdout():
    cmd = AtomicCmd(("echo", "Oil. Oil."), TEMP_OUT_STDOUT="dm")
    assert_equal(pformat(cmd), ("<Command = ['echo', 'Oil. Oil.']\n"
                                " STDOUT* = '${TEMP_DIR}/dm'\n"
                                " STDERR* = '${TEMP_DIR}/pipe_echo_%i.stderr'>") % (id(cmd),))


def test_pformat__atomiccmd__simple_with_temp_stdout__set_cwd():
    cmd = AtomicCmd(("echo", "Room service. Room service."),
                    TEMP_OUT_STDOUT="pv",
                    set_cwd=True)
    assert_equal(pformat(cmd), ("<Command = ['echo', 'Room service. Room service.']\n"
                                " STDOUT* = 'pv'\n"
                                " STDERR* = 'pipe_echo_%i.stderr'\n"
                                " CWD     = '${TEMP_DIR}'>") % (id(cmd),))


def test_pformat__atomiccmd__simple_with_stdout_pipe():
    cmd = AtomicCmd(("echo", "!"), OUT_STDOUT=AtomicCmd.PIPE)
    assert_equal(pformat(cmd), ("<Command = ['echo', '!']\n"
                                " STDOUT  = <PIPE>\n"
                                " STDERR* = '${TEMP_DIR}/pipe_echo_%i.stderr'>") % (id(cmd),))


###############################################################################
###############################################################################
# ParallelCmds

def test_pformat__sets__simple():
    def _do_test_pformat__sets__simple(cls, description):
        cmd_1 = AtomicCmd(("echo", "foo"), OUT_STDOUT=AtomicCmd.PIPE)
        cmd_2 = AtomicCmd("gzip", IN_STDIN=cmd_1)
        cmd = cls((cmd_1, cmd_2))
        assert_equal(pformat(cmd),
                     ("<{description}:\n"
                      "  - <00> Command = ['echo', 'foo']\n"
                      "         STDOUT  = <01>\n"
                      "         STDERR* = '${{TEMP_DIR}}/pipe_echo_{cmd_1_id}.stderr'\n"
                      "  - <01> Command = ['gzip']\n"
                      "         STDIN   = <00>\n"
                      "         STDOUT* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2_id}.stdout'\n"
                      "         STDERR* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2_id}.stderr'>")
                     .format(description=description,
                             cmd_1_id=id(cmd_1),
                             cmd_2_id=id(cmd_2)))
    yield _do_test_pformat__sets__simple, ParallelCmds, "Parallel commands"
    yield _do_test_pformat__sets__simple, SequentialCmds, "Sequential commands"


def test_pformat__sets__nested():
    cmd_1 = AtomicCmd(("echo", "foo"), OUT_STDOUT=AtomicCmd.PIPE)
    cmd_2 = AtomicCmd("gzip", IN_STDIN=cmd_1)
    cmd_3 = AtomicCmd("sha1sum")
    set_1 = ParallelCmds((cmd_1, cmd_2))
    set_2 = SequentialCmds((set_1, cmd_3))
    assert_equal(pformat(set_2),
                 ("<Sequential commands:\n"
                  "  - Parallel commands:\n"
                  "      - <00> Command = ['echo', 'foo']\n"
                  "             STDOUT  = <01>\n"
                  "             STDERR* = '${{TEMP_DIR}}/pipe_echo_{cmd_1_id}.stderr'\n"
                  "      - <01> Command = ['gzip']\n"
                  "             STDIN   = <00>\n"
                  "             STDOUT* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2_id}.stdout'\n"
                  "             STDERR* = '${{TEMP_DIR}}/pipe_gzip_{cmd_2_id}.stderr'\n"
                  "  - <02> Command = ['sha1sum']\n"
                  "         STDOUT* = '${{TEMP_DIR}}/pipe_sha1sum_{cmd_3_id}.stdout'\n"
                  "         STDERR* = '${{TEMP_DIR}}/pipe_sha1sum_{cmd_3_id}.stderr'>")
                 .format(cmd_1_id=id(cmd_1),
                         cmd_2_id=id(cmd_2),
                         cmd_3_id=id(cmd_3)))


###############################################################################
###############################################################################
# Bad input

def test_pformat__bad_input():
    @nose.tools.raises(TypeError)
    def _do_test_pformat__bad_input(value):
        pformat(value)
    yield _do_test_pformat__bad_input, 1
    yield _do_test_pformat__bad_input, {}
    yield _do_test_pformat__bad_input, ""


###############################################################################
###############################################################################
# _pformat_list

def test_pformat_list__empty():
    assert_equal(_pformat_list([]), "[]")


def test_pformat_list__single():
    assert_equal(_pformat_list([3]), "[3]")


def test_pformat_list__multiple():
    assert_equal(_pformat_list([3, 2, 1]), "[3, 2, 1]")


def test_pformat_list__wrapped():
    assert_equal(_pformat_list([3, 2, 1], width=1), "[3,\n 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=2), "[3,\n 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=3), "[3,\n 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=4), "[3,\n 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=5), "[3,\n 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=6), "[3, 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=7), "[3, 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=8), "[3, 2,\n 1]")
    assert_equal(_pformat_list([3, 2, 1], width=9), "[3, 2, 1]")
    assert_equal(_pformat_list([3, 2, 1], width=10), "[3, 2, 1]")
