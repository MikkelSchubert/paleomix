# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import string
import sys
from pathlib import Path
from typing import cast

from paleomix.common.argparse import ArgumentParser, Namespace


def parse_args(argv: list[str]) -> Namespace:
    parser = ArgumentParser("paleomix ngs:mask_fasta")
    parser.add_argument("fasta", type=Path)

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    table = list(range(256))
    for char in string.ascii_letters:
        if char not in "acgtnACGTN":
            table[ord(char)] = ord("N")

    table = bytes(table)

    buffer = sys.stdout.buffer
    with cast("Path", args.fasta).open("rb") as handle:
        for line in handle:
            if not line.startswith(b">"):
                line = line.translate(table)

            buffer.write(line)

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv[1:]))
