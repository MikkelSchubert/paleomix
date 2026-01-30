# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import os
import re
import subprocess
import sys

from paleomix.common import resources
from paleomix.common.argparse import ArgumentParser, Namespace


def parse_args(argv: list[str]) -> Namespace:
    parser = ArgumentParser(prog="paleomix :rscript")
    parser.add_argument("script", help="Name of Rscript bundled with PALEOMIX")
    parser.add_argument("args", nargs="*")

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    if not re.match(r"^[a-z]+/[a-z]+\.r$", args.script, flags=re.IGNORECASE):
        print(f"ERROR: Invalid resource {args.script!r}")
        return 1

    with resources.access(os.path.join("rscripts", args.script)) as path:
        result = subprocess.run(
            ("Rscript", path, *args.args),
            stdin=subprocess.DEVNULL,
            check=False,
        )

    return result.returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
