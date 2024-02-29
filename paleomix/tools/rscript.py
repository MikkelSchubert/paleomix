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
import os
import re
import subprocess
import sys

import paleomix.common.argparse as argparse
from paleomix.common import resources


def parse_args(argv):
    parser = argparse.ArgumentParser(prog="paleomix :rscript")
    parser.add_argument("script", help="Name of Rscript bundled with PALEOMIX")
    parser.add_argument("args", nargs="*")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    if not re.match(r"^[a-z]+/[a-z]+\.r$", args.script, flags=re.I):
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
