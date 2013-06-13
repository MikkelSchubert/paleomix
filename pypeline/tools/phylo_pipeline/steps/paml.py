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
from pypeline.node import Node, MetaNode, CommandNode, NodeError
from pypeline.nodes.sequences import CollectSequencesNode
from pypeline.nodes.mafft import MetaMAFFTNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicset import ParallelCmds

from pypeline.common.utilities import fragment, safe_coerce_to_tuple
from pypeline.common.formats.msa import read_msa
import pypeline.common.fileutils as fileutils


import common


class FastaToPAMLPhyNode(Node):
    def __init__(self, input_file, output_file, exclude_groups, dependencies = ()):
        self._input_file  = input_file
        self._output_file = output_file
        self._excluded = safe_coerce_to_tuple(exclude_groups)
        description  = "<FastaToPAMLPhy: '%s' -> '%s'>" % \
            (input_file, output_file)

        Node.__init__(self,
                      description  = description,
                      input_files  = [input_file],
                      output_files = [output_file],
                      dependencies = dependencies)


    def _run(self, _config, temp):
        msa = read_msa(self._input_file)
        for excluded_group in self._excluded:
            msa.pop(excluded_group)

        lines = []
        lines.append("  %i %i" % (len(msa), len(msa.itervalues().next())))
        for (name, seq) in sorted(msa.iteritems()):
            lines.append("")
            lines.append(name)

            for line in fragment(60, seq.upper()):
                lines.append(" ".join(fragment(3, line)))

        with open(fileutils.reroot_path(temp, self._output_file), "w") as output:
            output.write("\n".join(lines))


    def _teardown(self, _config, temp):
        source_file = fileutils.reroot_path(temp, self._output_file)
        output_file = self._output_file
        fileutils.move_file(source_file, output_file)


class CodemlNode(CommandNode):
    def __init__(self, control_file, sequence_file, trees_file, output_prefix, dependencies = ()):
        self._control_file  = control_file
        self._sequence_file = sequence_file
        self._trees_file    = trees_file
        self._output_prefix = output_prefix

        command = AtomicCmd(["codeml", "template.ctl"],
                            IN_CONTROL_FILE  = control_file,
                            IN_SEQUENCE_FILE = sequence_file,
                            IN_TREES_FILE    = trees_file,
                            TEMP_OUT_CTL     = "template.ctl",
                            TEMP_OUT_SEQS    = "template.seqs",
                            TEMP_OUT_TREES   = "template.trees",
                            TEMP_OUT_STDOUT  = "template.stdout",
                            TEMP_OUT_STDERR  = "template.stderr",
                            OUT_CODEML       = output_prefix + ".codeml",
                            TEMP_OUT_2NG_DN  = output_prefix + ".2NG.dN",
                            TEMP_OUT_2NG_DS  = output_prefix + ".2NG.dS",
                            TEMP_OUT_2NG_T   = output_prefix + ".2NG.t",
                            TEMP_OUT_4FOLD   = output_prefix + ".4fold.nuc",
                            TEMP_OUT_LNF     = output_prefix + ".lnf",
                            TEMP_OUT_RST     = output_prefix + ".rst",
                            TEMP_OUT_RST1    = output_prefix + ".rst1",
                            TEMP_OUT_RUB     = output_prefix + ".rub",
                            IN_STDIN         = "/dev/null", # Prevent promts from blocking
                            set_cwd          = True)

        CommandNode.__init__(self,
                             description  = "<CodemlNode: '%s' -> '%s.*'>" % (sequence_file, output_prefix),
                             command      = command,
                             dependencies = dependencies)

    def _setup(self, _config, temp):
        def _symlink(filename, dst):
            os.symlink(os.path.realpath(filename), os.path.join(temp, dst))
        _symlink(self._sequence_file, "template.seqs")
        _symlink(self._trees_file,    "template.trees")

        self._update_ctl_file(source        = self._control_file,
                              destination   = os.path.join(temp, "template.ctl"),
                              sequence_file = os.path.join(temp, "template.seqs"),
                              trees_file    = os.path.join(temp, "template.trees"),
                              output_prefix = self._output_prefix)

    def _run(self, config, temp):
        try:
            CommandNode._run(self, config, temp)
        except NodeError, error:
            # Allow failures due to low coverage
            with open(fileutils.reroot_path(temp, "template.stdout")) as handle:
                codeml = handle.read()
                if "sequences do not have any resolved nucleotides. Giving up." not in codeml:
                    raise error

            with open(fileutils.reroot_path(temp, self._output_prefix + ".codeml"), "a") as handle:
                handle.write("\nWARNING: No resolved nucleotides found, could not process gene.\n")

            import sys
            sys.stderr.write("WARNING: No resolved nucleotides in " + self._output_prefix + "\n")


    def _teardown(self, config, temp):
        prefix = os.path.basename(self._output_prefix)
        for filename in ("2NG.dN", "2NG.dS", "2NG.t", "4fold.nuc", "lnf", "rst", "rst1", "rub"):
            src_path = os.path.join(temp, filename)
            dst_path = os.path.join(temp, "%s.%s" % (prefix, filename))

            if os.path.exists(src_path):
                os.rename(src_path, dst_path)

        CommandNode._teardown(self, config, temp)

    @classmethod
    def _update_ctl_file(cls, source, destination, sequence_file, trees_file, output_prefix):
        output_file = output_prefix + ".codeml"
        with open(source) as handle:
            template = handle.read()

        # TODO: Check that number of replacements == 1
        # TODO: Do check before running everything!
        template, count = re.subn(r'(\s*seqfile\s*=).*',  r'\1 ' + os.path.basename(sequence_file), template)
        assert count == 1, count
        template, count = re.subn(r'(\s*treefile\s*=).*', r'\1 ' + os.path.basename(trees_file),    template)
        assert count == 1, count
        template, count = re.subn(r'(\s*outfile\s*=).*',  r'\1 ' + os.path.basename(output_file),   template)
        assert count == 1, count

        with open(destination, "w") as handle:
            handle.write(template)


def build_codeml_nodes(options, settings, interval, taxa, filtering, dependencies):
    in_postfix, out_postfix, afa_ext = "", "", ".afa"
    if any(filtering.itervalues()):
        in_postfix = out_postfix = ".filtered"
    if not settings["MSAlignment"]["Enabled"]:
        out_postfix = ".unaligned" + out_postfix
        afa_ext = ".fasta"

    paml        = settings["PAML"]
    sequences   = common.collect_sequences(options, interval, taxa)
    sequencedir = os.path.join(options.destination, "alignments", interval["Name"] + in_postfix)
    destination = os.path.join(options.destination, "paml", "codeml", interval["Name"] + out_postfix)


    # Build meta-node for sequence conversion to PHYLIP format accepted by codeml
    phylip_nodes = {}
    for sequence in sequences:
        input_file  = os.path.join(sequencedir, sequence + afa_ext)
        output_file = os.path.join(destination, sequence + ".phy")

        phylip_nodes[sequence] = FastaToPAMLPhyNode(input_file     = input_file,
                                                    output_file    = output_file,
                                                    exclude_groups = paml["codeml"]["ExcludeGroups"],
                                                    dependencies   = dependencies)

    phylip_meta = MetaNode(description  = "<FastaToPAMLPhyNodes: '%s/*.%s' -> '%s/*.phy'>" \
                           % (sequencedir, afa_ext, destination),
                           subnodes     = phylip_nodes.values(),
                           dependencies = dependencies)

    codeml_nodes = []
    for (ctl_name, ctl_file) in paml["codeml"]["Control Files"].iteritems():
        for (sequence, node) in phylip_nodes.iteritems():
            output_prefix = os.path.join(destination, sequence + ".%s" % (ctl_name,))

            codeml = CodemlNode(control_file  = ctl_file,
                                trees_file    = paml["codeml"]["Tree File"],
                                sequence_file = iter(node.output_files).next(),
                                output_prefix = output_prefix,
                                dependencies  = node)
            codeml_nodes.append(codeml)

    return MetaNode(description  = "<CodemlNodes>",
                    subnodes     = codeml_nodes,
                    dependencies = phylip_meta)



def chain_codeml(pipeline, options, makefiles):
    destination = options.destination # Move to makefile
    for makefile in makefiles:
        nodes     = []
        intervals = makefile["Project"]["Intervals"]
        filtering = makefile["Project"]["Filter Singletons"]
        taxa      = makefile["Project"]["Taxa"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        for interval in intervals.itervalues():
            nodes.append(build_codeml_nodes(options, makefile, interval, taxa, filtering, makefile["Nodes"]))
        makefile["Nodes"] = tuple(nodes)
    options.destination = destination



if __name__ == '__main__':
    import sys
    class config:
        temp_root = "./temp"
    
    node = CodemlNode(control_file  = sys.argv[1],
                      sequence_file = sys.argv[2],
                      trees_file    = sys.argv[3],
                      output_prefix = sys.argv[4])
    node.run(config)
