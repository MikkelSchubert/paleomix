# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import sys


def test_require_dev_mode() -> None:
    assert sys.flags.dev_mode
