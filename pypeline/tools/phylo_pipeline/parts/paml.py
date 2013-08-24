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

from pypeline.node import Node, MetaNode, CommandNode, NodeError
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import SequentialCmds


from pypeline.common.utilities import fragment, safe_coerce_to_tuple
from pypeline.common.formats.msa import MSA
import pypeline.common.fileutils as fileutils
import pypeline.tools.phylo_pipeline.parts.common as common


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
        msa = MSA.from_file(self._input_file)
        if self._excluded:
            msa = msa.exclude(self._excluded)

        lines = []
        lines.append("  %i %i" % (len(msa), msa.seqlen()))
        for record in sorted(msa):
            lines.append("")
            lines.append(record.name)

            for line in fragment(60, record.sequence.upper()):
                lines.append(" ".join(fragment(3, line)))

        with open(fileutils.reroot_path(temp, self._output_file), "w") as output:
            output.write("\n".join(lines))


    def _teardown(self, _config, temp):
        source_file = fileutils.reroot_path(temp, self._output_file)
        output_file = self._output_file
        fileutils.move_file(source_file, output_file)


class CodemlNode(CommandNode):
    def __init__(self, control_file, sequence_file, trees_file, output_tar, dependencies = ()):
        self._control_file  = control_file
        self._sequence_file = sequence_file
        self._trees_file    = trees_file

        paml_cmd = AtomicCmd(["codeml", "template.ctl"],
                             IN_CONTROL_FILE  = control_file,
                             IN_SEQUENCE_FILE = sequence_file,
                             IN_TREES_FILE    = trees_file,
                             TEMP_OUT_CTL     = "template.ctl",
                             TEMP_OUT_SEQS    = "template.seqs",
                             TEMP_OUT_TREES   = "template.trees",
                             TEMP_OUT_STDOUT  = "template.stdout",
                             TEMP_OUT_STDERR  = "template.stderr",
                             IN_STDIN         = "/dev/null", # Prevent promts from blocking
                             set_cwd          = True,
                             **CodemlNode._get_codeml_files("TEMP_OUT_CODEML"))

        tar_pairs = CodemlNode._get_codeml_files("TEMP_IN_CODEML")
        tar_files = ["%%(%s)s" % (key,) for key in tar_pairs]
        tar_cmd  = AtomicCmd(["tar", "cvzf", "%(OUT_FILE)s"] + tar_files,
                             OUT_FILE = output_tar,
                             set_cwd  = True,
                             **tar_pairs)

        CommandNode.__init__(self,
                             description  = "<CodemlNode: %r -> %r>" % (sequence_file, output_tar),
                             command      = SequentialCmds([paml_cmd, tar_cmd]),
                             dependencies = dependencies)


    def _setup(self, _config, temp):
        def _symlink(filename, dst):
            os.symlink(os.path.realpath(filename), os.path.join(temp, dst))
        _symlink(self._sequence_file, "template.seqs")
        _symlink(self._trees_file,    "template.trees")

        self._update_ctl_file(source      = self._control_file,
                              destination = os.path.join(temp, "template.ctl"))


    def _run(self, config, temp):
        try:
            CommandNode._run(self, config, temp)
        except NodeError, error:
            if self._command.join() == [1, None]:
                with open(fileutils.reroot_path(temp, "template.stdout")) as handle:
                    lines = handle.readlines()
                if lines and ("Giving up." in lines[-1]):
                    error = NodeError("%s\n\n%s" % (error, lines[-1]))
            raise error


    @classmethod
    def _update_ctl_file(cls, source, destination):
        with open(source) as handle:
            template = handle.read()

        # TODO: Do check before running everything!
        template, count = re.subn(r'(\bseqfile\s*=).*',  r'\1 template.seqs', template)
        assert count == 1, count
        template, count = re.subn(r'(\btreefile\s*=).*', r'\1 template.trees', template)
        assert count == 1, count
        template, count = re.subn(r'(\boutfile\s*=).*',  r'\1 mlc',       template)
        assert count == 1, count

        with open(destination, "w") as handle:
            handle.write(template)

    @classmethod
    def _get_codeml_files(cls, key_type):
        results      = {}
        codeml_files = ["mlc", "2NG.dN", "2NG.dS", "2NG.t", "4fold.nuc", "lnf", "rst", "rst1", "rub"]
        for filename in codeml_files:
            key = "%s_%s" % (key_type, filename.upper().replace(".", "_"))
            results[key] = filename
        return results


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
    for (ctl_name, ctl_files) in paml["codeml"].iteritems():
        if ctl_name in ("ExcludeGroups",):
            continue

        for (sequence, node) in phylip_nodes.iteritems():
            output_tar = os.path.join(destination, "%s.%s.tar.gz" % (sequence, ctl_name))

            codeml = CodemlNode(control_file  = ctl_files["Control File"],
                                trees_file    = ctl_files["Tree File"],
                                sequence_file = iter(node.output_files).next(),
                                output_tar    = output_tar,
                                dependencies  = node)
            codeml_nodes.append(codeml)

    return MetaNode(description  = "<CodemlNodes>",
                    subnodes     = codeml_nodes,
                    dependencies = phylip_meta)



def chain_codeml(_pipeline, options, makefiles):
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
