#!/usr/bin/python
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
from paleomix.node import \
     Node
from paleomix.common.fileutils import \
     copy_file, \
     reroot_path
from paleomix.common.utilities import \
     safe_coerce_to_tuple


class CopyOutputFilesNode(Node):
    """Copies the output-files of one or more nodes to a specified folder."""

    def __init__(self, description, destination, source_nodes):
        source_nodes = safe_coerce_to_tuple(source_nodes)

        input_files  = []
        for source_node in source_nodes:
            input_files.extend(source_node.output_files)

        output_files = [reroot_path(destination, fpath) for fpath in input_files]
        self._files  = zip(input_files, output_files)

        Node.__init__(self,
                      description  = "<Copy %s output to %r>" % (description, destination),
                      input_files  = input_files,
                      output_files = output_files,
                      dependencies = source_nodes)


    def _run(self, _config, _temp):
        for (src_file, dst_file) in self._files:
            copy_file(src_file, dst_file)
