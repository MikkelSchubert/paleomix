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

import signal

from paleomix.common.procs import (
    RegisteredPopen,
    running_processes,
    terminate_all_processes,
    terminate_processes,
)


# Test that the internal list of processes is kept clean of joined processes objects
def test_atomiccmd__cleanup_proc__commit() -> None:
    proc = RegisteredPopen("true")
    assert proc in running_processes()
    assert proc.wait() == 0
    assert proc not in running_processes()


def test_atomiccmd__cleanup_proc__gc() -> None:
    proc = RegisteredPopen("true")
    assert proc in running_processes()

    pid = proc.pid
    del proc

    try:
        assert any(proc.pid == pid for proc in running_processes())
    finally:
        terminate_all_processes()


def test_atomiccmd__terminate_all_processes() -> None:
    proc = RegisteredPopen(("sleep", "10"))
    assert proc in running_processes()
    terminate_all_processes()
    assert proc not in running_processes()
    assert proc.poll() == -signal.SIGTERM


def test_atomiccmd__terminate_processes() -> None:
    proc_1 = RegisteredPopen(("sleep", "10"))
    proc_2 = RegisteredPopen(("sleep", "10"))

    assert proc_1 in running_processes()
    assert proc_2 in running_processes()

    terminate_processes([proc_1])

    assert proc_1 not in running_processes()
    assert proc_2 in running_processes()

    assert proc_1.poll() == -signal.SIGTERM
    assert proc_2.poll() is None

    terminate_all_processes()
