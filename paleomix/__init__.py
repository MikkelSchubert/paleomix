#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
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
#

__version_info__ = (1, 2, 6)
__version__ = '%i.%i.%i' % __version_info__


def run(command=None):
    """Main entry-point for setuptools"""
    import sys
    import paleomix.main

    argv = []
    if command is not None:
        argv.append(command)
    argv.extend(sys.argv[1:])

    return paleomix.main.main(argv)


def run_bam_pipeline():
    """Legacy entry-point for setuptools"""
    return run("bam_pipeline")


def run_gtf_to_bed():
    """Legacy entry-point for setuptools"""
    return run("gtf_to_bed")


def run_phylo_pipeline():
    """Legacy entry-point for setuptools"""
    return run("phylo_pipeline")


def run_rmdup_collapsed():
    """Legacy entry-point for setuptools"""
    return run("rmdup_collapsed")


def run_trim_pipeline():
    """Legacy entry-point for setuptools"""
    return run("trim_pipeline")
