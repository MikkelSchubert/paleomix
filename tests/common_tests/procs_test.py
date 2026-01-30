# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
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
