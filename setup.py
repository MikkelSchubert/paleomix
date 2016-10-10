#!/usr/bin/env python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import codecs
import os
import sys

from setuptools import setup, find_packages


if (sys.version_info[0] != 2) or (sys.version_info[1] != 7):
    sys.stderr.write("ERROR: Python version 2.7.x required!\n")
    sys.stderr.write("       Current version is v%s\n"
                     % (sys.version.replace("\n", " "),))
    sys.exit(1)


def _get_version():
    """Retrieve version from current install directory."""
    env = {}
    with open(os.path.join("paleomix", "__init__.py")) as handle:
        exec(handle.read(), env)

    return env["__version__"]


def _get_readme():
    """Retrieves contents of README.rst, forcing UTF-8 encoding."""
    with codecs.open("README.rst", encoding="utf-8") as handle:
        return handle.read()


setup(
    name='paleomix',
    version=_get_version(),

    description='Bioinformatics pipelines for HTS data',
    long_description=_get_readme(),

    url='https://github.com/MikkelSchubert/paleomix',

    author='Mikkel Schubert',
    author_email='MSchubert@snm.ku.dk',

    license='MIT',

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2 :: Only',
        'Programming Language :: Python :: 2.7',
    ],

    keywords='pipeline bioinformatics hts phylogeny bam',

    packages=find_packages(exclude=['misc', 'tests']),

    install_requires=['pysam>=0.8.3',
                      'setproctitle>=1.1.0'],

    # Dependencies set in setup_requires to allow use of 'setup.py nosetests'
    setup_requires=['nose>=1.3.0',
                    'flexmock>=0.9.7',
                    'coverage>=4.0.0'],

    test_suite='nose.collector',

    entry_points={
        'console_scripts': [
            'paleomix=paleomix:run',

            # Aliases used in previous publications
            'bam_pipeline=paleomix:run_bam_pipeline',
            'bam_rmdup_collapsed=paleomix:run_rmdup_collapsed',
            'conv_gtf_to_bed=paleomix:run_gtf_to_bed',
            'phylo_pipeline=paleomix:run_phylo_pipeline',
            'trim_pipeline=paleomix:run_trim_pipeline',
        ],
    },

    zip_safe=False,

    include_package_data=True,
)
