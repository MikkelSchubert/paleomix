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
import re
import glob

from paleomix.common.makefile import \
    MakefileError


def is_paired_end(template):
    """Returns true if a template contains a Pair component."""
    return (template.format(Pair=1) != template)


def collect_files(path, template):
    """

    """
    if is_paired_end(template):
        if _has_glob_magic(template):
            result = {"PE_1": _sorted_glob(template.format(Pair=1)),
                      "PE_2": _sorted_glob(template.format(Pair=2))}

            if not (result["PE_1"] or result["PE_2"]):
                _raise_missing_files("paired-end", path, template)
            elif len(result["PE_1"]) != len(result["PE_2"]):
                raise MakefileError("Unequal number of mate 1 and mate 2 "
                                    "files found at path %r; found %i mate 1 "
                                    "files, and %i mate 2 files; specified in "
                                    "makefile at %r. Please verify that the "
                                    "path is correct, and update the makefile!"
                                    % (template,
                                       len(result["PE_1"]),
                                       len(result["PE_2"]),
                                       " :: ".join(path)))
        else:
            result = {"PE_1": [template.format(Pair=1)],
                      "PE_2": [template.format(Pair=2)]}
    elif _has_glob_magic(template):
        result = {"SE": _sorted_glob(template)}
        if not result["SE"]:
            _raise_missing_files("single-end", path, template)
    else:
        result = {"SE": [template]}

    return result


def _has_glob_magic(filename):
    """Returns true if a path contains wildcards used by glob / fnmatch."""
    return _GLOB_MAGIC.search(filename) is not None


def _sorted_glob(tmpl):
    return list(sorted(glob.iglob(tmpl)))


def _raise_missing_files(description, path, template):
    raise MakefileError("No files found for %s reads using path %r; "
                        "specified in makefile at %r. Please verify that the "
                        "path is correct, and update the makefile!"
                        % (description, template, " :: ".join(path)))


_GLOB_MAGIC = re.compile('[*?[]')
