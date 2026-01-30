# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import os

from paleomix.common import versions
from paleomix.tools.factory import rscript_command


def requirement(
    module: str,
    specifiers: str = "",
) -> versions.Requirement:
    return versions.Requirement(
        call=rscript_command((os.path.join("common", "requires.r"), module)),
        # d0fd3ea6 is a magic value printed by the requires.r script. This is used to
        # ensure that we can differentiate between the version and any other output
        regexp=r"d0fd3ea6: (\d+)\.(\d+)(?:\.(\d+))?",
        specifiers=specifiers,
        name=f"R module: {module}",
    )
