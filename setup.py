# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
#!/usr/bin/python
import os
import re
import sys
import subprocess

from distutils.core import setup
from distutils.command.build_py import build_py


if (sys.version_info.major != 2) or (sys.version_info.minor != 7):
    sys.stderr.write("ERROR: Python version 2.7.x required!\n")
    sys.stderr.write("       Current version is v%s\n" % (sys.version.replace("\n", " "),))
    sys.exit(1)


_VERSION = None
def get_version():
    global _VERSION
    if _VERSION is None:
        command = ("git", "describe", "--always", "--tags", "--dirty")
        try:
            _VERSION = subprocess.check_output(command).strip()
        except (subprocess.CalledProcessError, OSError), error:
            raise SystemExit(("Could not determine pypeline version:\n"
                              "  Command = %s\n"
                              "  Error   = %s") % (" ".join(command), error))
    return _VERSION


def locate_packages():
    packages = ['pypeline']
    for (dirpath, dirnames, _) in os.walk(packages[0]):
        for dirname in dirnames:
            package = os.path.join(dirpath, dirname).replace(os.sep, ".")
            packages.append(package)
    return packages


def locate_scripts():
    scripts = []
    for filename in os.listdir("bin"):
        if re.match(r"^[0-9a-z_]+", filename):
            script = os.path.join("bin", filename)
            if os.path.isfile(script) or os.path.islink(script):
                scripts.append(script)
    return scripts


class setup_version(build_py):
    def run(self):
        version = get_version()
        with open(os.path.join("pypeline", "version.py"), "w") as handle:
            handle.write("#!/usr/bin/env python\n")
            handle.write("__version__ = %r\n" % (version.strip(),))
        build_py.run(self)


setup(name         = 'Pypeline',
      version      = get_version(),
      description  = '(Framework for) Bioinformatics pipelines',
      author       = 'Mikkel Schubert',
      author_email = 'MSchubert@snm.ku.dk',
      url          = 'https://github.com/MikkelSchubert/pypeline',
      requires     = ['pysam (>=0.7.4)'],
      packages     = locate_packages(),
      scripts      = locate_scripts(),
      cmdclass     = {"build_py" : setup_version},
    )
