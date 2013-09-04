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
from pypeline.common.utilities import safe_coerce_to_tuple
from pypeline.node import MetaNode

class Target:
    def __init__(self, config, prefixes, name):
        self.name     = name
        self.prefixes = safe_coerce_to_tuple(prefixes)

        self._nodes_extras    = {}
        self._nodes_alignment = MetaNode(description  = "Alignments:",
                                          dependencies = [prefix.node for prefix in self.prefixes])

        self._setup_mapdamage_nodes()


    def add_extra_nodes(self, name, nodes):
        if name in self._nodes_extras:
            raise RuntimeError("Cannot add extra nodes '%s', task with that name exists!" % name)

        self._nodes_extras[name] = MetaNode(description  = name + ":",
                                            dependencies = nodes)


    @property
    def node(self):
        extras = MetaNode(description  = "Additional tasks:",
                          dependencies = self._nodes_extras.values())

        return MetaNode(description    = "Target: %s" % self.name,
                        dependencies   = (self._nodes_alignment, extras))


    def _setup_mapdamage_nodes(self):
        mapdamage_nodes = []
        for prefix in self.prefixes:
            prefix_nodes = []
            for sample in prefix.samples:
                for library in sample.libraries:
                    if library.mapdamage:
                        prefix_nodes.append(library.mapdamage)

            if any(prefix_nodes):
                node = MetaNode(description = prefix.name,
                                subnodes    = prefix_nodes)
                mapdamage_nodes.append(node)

        if mapdamage_nodes:
            self.add_extra_nodes("mapDamage", mapdamage_nodes)
