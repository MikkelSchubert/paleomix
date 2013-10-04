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
import logilab.astng
from logilab.astng import MANAGER


_DECORATOR = 'pypeline.atomiccmd.builder.create_customizable_cli_parameters'


def _disable_infer(_self, *_args, **_kwargs):
    raise logilab.astng.InferenceError()


def atomiccmd_builder_transform(module):
    if isinstance(module, logilab.astng.Function):
        if _DECORATOR in module.decoratornames():
            if module.type == 'method':
                module.type = 'classmethod'

            # Crude workaround to spurious errors:
            # Prevent pylint from attempting to infer the return type
            module.infer_call_result = _disable_infer
            return

    for child in module.get_children():
        atomiccmd_builder_transform(child)


def register(*_args, **_kwargs):
    MANAGER.register_transformer(atomiccmd_builder_transform)
