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
import re

from pypeline import Pypeline
from pypeline.node import Node, MetaNode, CommandNode
from pypeline.nodes.sequences import CollectSequencesNode, \
    FilterSingletonsMetaNode
from pypeline.nodes.mafft import MetaMAFFTNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds

from pypeline.common.utilities import fragment
from pypeline.common.formats.msa import read_msa
import pypeline.common.fileutils as fileutils


import common


class FastaToPAMLPhyNode(Node):
    def __init__(self, input_file, output_file, dependencies = ()):
        description  = "<FastaToPAMLPhy: '%s' -> '%s'>" % \
            (input_file, output_file)
            
        Node.__init__(self, 
                      description  = description,
                      input_files  = [input_file],
                      output_files = [output_file],
                      dependencies = dependencies)


    def _run(self, _config, temp):
        msa = read_msa(self.input_files[0])
        
        lines = []
        lines.append("  %i %i" % (len(msa), len(msa.itervalues().next())))
        for (name, seq) in sorted(msa.iteritems()):
            lines.append("")
            lines.append(name)

            for line in fragment(60, seq.upper()):
                lines.append(" ".join(fragment(3, line)))

        with open(fileutils.reroot_path(temp, self.output_files[0]), "w") as output:
            output.write("\n".join(lines))


    def _teardown(self, _config, temp):
        output_file = self.output_files[0]
        fileutils.move_file(fileutils.reroot_path(temp, output_file), output_file)

        
class CodemlNode(CommandNode):
    def __init__(self, control_file, sequence_file, trees_file, output_file, dependencies = ()):
        self._control_file  = control_file
        self._sequence_file = sequence_file
        self._trees_file    = trees_file
        self._output_file   = output_file
        
        command = AtomicCmd(["codeml", "template.ctl"],
                            IN_SEQUENCE_FILE = sequence_file,
                            IN_TREES_FILE    = trees_file,
                            OUT_RESULTS      = output_file,
                            IN_STDIN         = "/dev/null") # Prevent promts from blocking

        CommandNode.__init__(self,
                             description  = "<CodemlNode: '%s' -> '%s'>" % (sequence_file, output_file),
                             command      = command,
                             dependencies = dependencies)

    def _setup(self, _config, temp):
        self._update_ctl_file(source        = self._control_file,
                              destination   = os.path.join(temp, "template.ctl"),
                              sequence_file = self._sequence_file,
                              trees_file    = self._trees_file,
                              output_file   = self._output_file)

    def _run(self, config, temp):
        temp = os.path.realpath(temp)
        oldwd = os.getcwd()
        os.chdir(temp)
        try:
            CommandNode._run(self, config, temp)
        finally:
            os.chdir(oldwd)

    def _teardown(self, config, temp):
        # Remove everything but the output file
        # Given the large number of temporary files, this is a lot easier than specifying them all ...
        output_file = os.path.basename(self._output_file)
        for filename in os.listdir(temp):
            if (filename != output_file) and not filename.startswith("pipe_"):
                os.remove(os.path.join(temp, filename))

        CommandNode._teardown(self, config, temp)
            
    @classmethod 
    def _update_ctl_file(cls, source, destination, sequence_file, trees_file, output_file):
        with open(source) as handle:
            template = handle.read()

        # TODO: Check that number of replacements == 1
        # TODO: Do check before running everything!
        template, count = re.subn(r'(\s*seqfile\s*=).*',  r'\1 ' + os.path.realpath(sequence_file), template)
        template, count = re.subn(r'(\s*treefile\s*=).*', r'\1 ' + os.path.realpath(trees_file),    template)
        template, count = re.subn(r'(\s*outfile\s*=).*',  r'\1 ' + os.path.basename(output_file),   template)

        with open(destination, "w") as handle:
            handle.write(template)

    

def build_codeml_nodes(options, settings, interval, taxa, filtering, dependencies):
    postfix = ""
    if any(filtering.itervalues()):
        postfix = ".filtered"

    sequences   = common.collect_sequences(options, interval, taxa)
    sequencedir = os.path.join(options.destination, "alignments", interval["Name"] + postfix)
    destination = os.path.join(options.destination, "paml", "codeml")

    # Build meta-node for sequence conversion to PHYLIP format accepted by codeml
    phylip_nodes = {}
    for sequence in sequences:
        input_file  = os.path.join(sequencedir, sequence + ".afa")
        output_file = os.path.join(destination, sequence + ".phy")
        
        phylip_nodes[sequence] = FastaToPAMLPhyNode(input_file, output_file, dependencies)

    return MetaNode(description  = "<FastaToPAMLPhyNodes: '%s/*.afa' -> '%s/*.phy'>" % (sequencedir, destination),
                    subnodes     = phylip_nodes.values(),
                    dependencies = dependencies)
                    
    


def chain_codeml(pipeline, options, makefiles):
    destination = options.destination # Move to makefile
    for makefile in makefiles:
        nodes     = []
        settings  = makefile["MSAlignment"]
        intervals = makefile["Project"]["Intervals"]
        filtering = makefile["Project"]["Filter Singletons"]
        taxa      = makefile["Project"]["Taxa"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        for interval in intervals.itervalues():
            nodes.append(build_codeml_nodes(options, settings, interval, taxa, filtering, makefile["Nodes"]))
        makefile["Nodes"] = tuple(nodes)
    options.destination = destination



if __name__ == '__main__':
    import sys
    class config:
        temp_root = "./temp"
    
    node = CodemlNode(control_file  = sys.argv[1],
                      sequence_file = sys.argv[2], 
                      trees_file    = sys.argv[3],
                      output_file   = sys.argv[4])
    node.run(config)
