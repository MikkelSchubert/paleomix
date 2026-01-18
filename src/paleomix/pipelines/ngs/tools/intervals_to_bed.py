# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import sys
from pathlib import Path
from typing import cast

from paleomix.common.argparse import ArgumentParser, Namespace


def parse_args(argv: list[str]) -> Namespace:
    parser = ArgumentParser("paleomix ngs:intervals_to_bed")
    parser.add_argument("intervals", type=Path)

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    buffer = sys.stdout.buffer
    with cast("Path", args.intervals).open("rb") as handle:
        for line in handle:
            if not line.startswith(b"@"):
                chrom, start, end, *_ = line.split(b"\t", 3)
                buffer.write(b"%s\t%s\t%s\n" % (chrom, start, end))

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv[1:]))
