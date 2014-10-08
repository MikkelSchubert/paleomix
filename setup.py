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
import shutil

from distutils.core import \
    setup
from distutils.command.install_scripts import \
    install_scripts as DistutilsInstallScripts


if (sys.version_info.major != 2) or (sys.version_info.minor != 7):
    sys.stderr.write("ERROR: Python version 2.7.x required!\n")
    sys.stderr.write("       Current version is v%s\n"
                     % (sys.version.replace("\n", " "),))
    sys.exit(1)


def get_version():
    import pypeline
    return pypeline.__version__


def locate_packages():
    packages = ['pypeline']
    for (dirpath, dirnames, _) in os.walk(packages[0]):
        for dirname in dirnames:
            package = os.path.join(dirpath, dirname).replace(os.sep, ".")
            packages.append(package)
    return packages


def locate_scripts(is_link=False):
    scripts = []
    for filename in os.listdir("bin"):
        if re.match(r"^[0-9a-z_]+", filename):
            script = os.path.join("bin", filename)
            if os.path.islink(script) == is_link:
                scripts.append(script)
    return scripts


class InstallLinks(DistutilsInstallScripts):
    def run(self):
        DistutilsInstallScripts.run(self)

        # Ensure that symbolic links are preserved as symbolic links
        links = [os.path.basename(fpath)
                 for fpath in locate_scripts(is_link=True)]

        for bin_fpath in self.get_outputs():
            if os.path.basename(bin_fpath) == 'paleomix':
                fpath = os.path.dirname(bin_fpath)

                for link in links:
                    link_fpath = os.path.join(fpath, link)
                    print "linking", link_fpath, "-> 'paleomix'"
                    if os.path.exists(link_fpath):
                        os.remove(link_fpath)

                    os.symlink('paleomix', link_fpath)
                    shutil.copymode(bin_fpath, link_fpath)


setup(name='PALEOMIX Pipeline',
      version=get_version(),
      description='(Framework for) Bioinformatics pipelines',
      author='Mikkel Schubert',
      author_email='MSchubert@snm.ku.dk',
      url='https://github.com/MikkelSchubert/paleomix',
      requires=['pysam (>=0.7.5)'],
      packages=locate_packages(),
      scripts=locate_scripts(),
      cmdclass={'install_scripts': InstallLinks})
