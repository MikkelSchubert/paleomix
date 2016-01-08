#!/usr/bin/python
#
# Copyright (c) 2016 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import paleomix.resources

import paleomix.common.versions as versions

from paleomix.resources import rscript


def requirement(module, checks=versions.Any(), cache={}):
    key = (module, checks)
    result = cache.get(key)
    if result is None:
        filename = rscript("common", "requires.r")
        result = versions.Requirement(call=("Rscript", filename, module),
                                      search="d0fd3ea6: (\d+)\.(\d+)(?:\.(\d+))?",
                                      checks=checks,
                                      name="R module: {}".format(module))

        cache[(module, checks)] = result

    return result
