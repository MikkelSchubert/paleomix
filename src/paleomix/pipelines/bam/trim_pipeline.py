# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import paleomix.pipelines.bam.main as bam_main


def main(argv: list[str]) -> int:
    """Wrapper to invoke the trimming pipeline; used by paleomix.main."""
    return bam_main.main(argv, pipeline="trim")
