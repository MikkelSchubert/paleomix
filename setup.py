#
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
import sys

from setuptools import find_packages, setup

if sys.version_info < (3, 7):  # noqa: UP036
    print(
        "ERROR: PALEOMIX requires at least Python 3.7, but setup.py was run using "
        "Python v{}.{}".format(*sys.version_info),
        file=sys.stderr,
    )
    sys.exit(1)


def _get_version() -> str:
    """Retrieve version from current install directory."""
    with open(os.path.join("paleomix", "__init__.py")) as handle:
        for line in handle:
            match = re.match(r'__version__ = "(.*)"', line)
            if match is not None:
                return match.group(1)

    raise AssertionError("Could not determine PALEOMIX version")


def _get_readme() -> str:
    """Retrieves contents of README.rst, forcing UTF-8 encoding."""
    with open("README.rst", encoding="utf-8") as handle:
        return handle.read()


setup(
    name="paleomix",
    version=_get_version(),
    description="Bioinformatics pipelines for HTS data",
    long_description=_get_readme(),
    long_description_content_type="text/x-rst",
    url="https://github.com/MikkelSchubert/paleomix",
    author="Mikkel Schubert",
    author_email="MikkelSch@gmail.com",
    license="MIT",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="pipeline bioinformatics hts bam",
    packages=find_packages(exclude=["misc", "tests"]),
    install_requires=[
        "coloredlogs>=10.0",
        "configargparse>=0.13.0",
        "humanfriendly>=4.7",
        "packaging>=19.0",
        "pysam>=0.10.0",
        "ruamel.yaml>=0.16.0",
        "setproctitle>=1.1.0",
        "typing_extensions>=4.0",
    ],
    entry_points={"console_scripts": ["paleomix=paleomix.__main__:main"]},
    include_package_data=True,
)
