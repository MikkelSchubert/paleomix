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

from pypeline.common.formats.newick import \
     Newick
from pypeline.common.utilities import \
    safe_coerce_to_tuple
from pypeline.common.fileutils import \
    describe_files, \
    move_file
from pypeline.node import \
    Node


class NewickSupportNode(Node):
    def __init__(self, main_tree_files, support_tree_files, output_file, dependencies = ()):
        self._output_file        = output_file
        self._main_tree_files    = safe_coerce_to_tuple(main_tree_files)
        self._support_tree_files = safe_coerce_to_tuple(support_tree_files)
        input_files = self._main_tree_files + self._support_tree_files

        description  = "<NewickSupport: %s>" % \
          (describe_files(main_tree_files),)

        Node.__init__(self,
                      description  = description,
                      input_files  = input_files,
                      output_files = output_file,
                      dependencies = dependencies)

    def _run(self, _config, temp):
        main_trees    = self._read_trees(self._main_tree_files)
        support_trees = self._read_trees(self._support_tree_files)

        lines = []
        for main_tree in main_trees:
            rooted_tree    = main_tree.reroot_on_midpoint()
            supported_tree = rooted_tree.add_support(support_trees)
            lines.append(str(supported_tree))
        lines = "\n".join(lines) + "\n"

        temp_output_file = os.path.join(temp, os.path.basename(self._output_file))
        with open(temp_output_file, "w") as handle:
            handle.write(lines)

        move_file(temp_output_file, self._output_file)


    @classmethod
    def _read_trees(cls, filenames):
        trees = []
        for filename in filenames:
            with open(filename) as handle:
                for line in handle:
                    trees.append(Newick.from_string(line))
        return trees
