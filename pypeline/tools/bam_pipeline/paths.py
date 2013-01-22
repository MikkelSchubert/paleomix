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
#
import os
import glob
import itertools
import collections


def is_paired_end(template):
    return (template.format(Pair = 1) != template)


def collect_files(template):
    """ """
    if not is_paired_end(template):
        return {"SE" : list(sorted(glob.glob(template)))}
    
    files = {}
    for (ii, name) in enumerate(("PE_1", "PE_2"), start = 1):
        files[name] = list(sorted(glob.glob(template.format(Pair = ii))))
    return files


def reference_sequence(bwa_prefix):
    extensions = (".fa", ".fasta")
    if any(bwa_prefix.endswith(ext) for ext in extensions):
        return bwa_prefix
    
    for ext in extensions:
        if os.path.exists(bwa_prefix + ext):
            return bwa_prefix + ext
