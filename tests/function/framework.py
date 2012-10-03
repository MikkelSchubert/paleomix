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
import sys
import time
import traceback

import pypeline.ui as ui
from pypeline.node import Node
from pypeline.pipeline import Pypeline


_epoch = time.time()
class Config:
    temp_root    = "tests/runs/%i/temp/" % _epoch
    destination  = "tests/runs/%i/output/" % _epoch
    dependencies = Node(description = "DummyDependency")


def run(*functions):
    ui.print_info("Running functional Node tests:")
    ui.print_info("  - Destination = '%s' ..." % Config.destination)
    ui.print_info("  - Temp root = '%s' ..." % Config.temp_root)
    ui.print_info()

    errors = False
    pipeline = Pypeline(Config)
    for func in functions:
        try: 
            ui.print_info("Adding node '%s' ..." % func.__name__)
            node = func(Config)
            if list(node.dependencies) != [Config.dependencies]:
                raise RuntimeError("Node did not pass dependencies")

            pipeline.add_nodes(node)
        except Exception, e:
            ui.print_err(traceback.format_exc())
            errors = True

    if not pipeline.run(dry_run = True, collapse = False):
        errors = True
    elif not pipeline.run(collapse = False):
        errors = True

            
    sys.exit(1 if errors else 0)
