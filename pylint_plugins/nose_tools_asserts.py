#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is herby granted, free of charge, to any person obtaining a copy
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
from logilab.astng import MANAGER, scoped_nodes

_ASSERTS = ['equal', 'not_equal',
            'in', 'not_in',
            'is', 'is_not',
            'is_none', 'is_not_none',
            'is_instance', 'is_not_instance',
            'true', 'false',
            'raises', 'raises_regexp']

def nose_tools_transform(module):
    if module.name == 'nose.tools':
        for assert_func in _ASSERTS:
            assert_func = "assert_" + assert_func
            module.locals[assert_func] = [scoped_nodes.Class(assert_func, None)]

def register(*_args, **_kwargs):
    MANAGER.register_transformer(nose_tools_transform)
