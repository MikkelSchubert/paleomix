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

from paleomix.node import CommandNode, NodeError
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import SequentialCmds
from paleomix.common.formats.fasta import FASTA
from paleomix.common.utilities import safe_coerce_to_frozenset

import paleomix.common.fileutils as fileutils


class CodemlNode(CommandNode):
    def __init__(self, control_file, sequence_file, trees_file, output_tar,
                 exclude_groups=(), dependencies=()):
        self._exclude_groups = safe_coerce_to_frozenset(exclude_groups)
        self._control_file = control_file
        self._sequence_file = sequence_file
        self._trees_file = trees_file

        paml_cmd = AtomicCmd(["codeml", "template.ctl"],
                             IN_CONTROL_FILE  = control_file,
                             IN_SEQUENCE_FILE = sequence_file,
                             IN_TREES_FILE    = trees_file,
                             TEMP_OUT_CTL     = "template.ctl",
                             TEMP_OUT_SEQS    = "template.seqs",
                             TEMP_OUT_TREES   = "template.trees",
                             TEMP_OUT_STDOUT  = "template.stdout",
                             TEMP_OUT_STDERR  = "template.stderr",
                             TEMP_OUT_4FOLD   = "4fold.nuc",
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
        self._update_ctl_file(source      = self._control_file,
                              destination = os.path.join(temp, "template.ctl"))

        os.symlink(os.path.abspath(self._trees_file), os.path.join(temp, "template.trees"))
        with open(os.path.join(temp, "template.seqs"), "w") as handle:
            for record in FASTA.from_file(self._sequence_file):
                if record.name not in self._exclude_groups:
                    name = record.name
                    sequence = record.sequence.upper()
                    handle.write("%s\n" % (FASTA(name, None, sequence),))

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
        results = {}
        codeml_files = ["mlc", "2NG.dN", "2NG.dS", "2NG.t",
                        "lnf", "rst", "rst1", "rub"]

        for filename in codeml_files:
            key = "%s_%s" % (key_type, filename.upper().replace(".", "_"))
            results[key] = filename
        return results


def build_codeml_nodes(options, settings, regions, filtering, dependencies):
    in_postfix, out_postfix, afa_ext = "", "", ".afa"
    if any(filtering.itervalues()):
        in_postfix = out_postfix = ".filtered"
    if not settings["MultipleSequenceAlignment"][regions["Name"]]["Enabled"]:
        out_postfix = ".unaligned" + out_postfix
        afa_ext = ".fasta"

    codeml = settings["PAML"]["codeml"]
    subset_key = codeml["SubsetRegions"].get(regions["Name"])
    sequences = regions["Sequences"][subset_key]
    sequencedir = os.path.join(options.destination, "alignments", regions["Name"] + in_postfix)
    destination = os.path.join(options.destination, "paml", "codeml", regions["Name"] + out_postfix)

    fasta_files = {}
    for sequence in sequences:
        fasta_files[sequence] = os.path.join(sequencedir, sequence + afa_ext)

    codeml_nodes = []
    for (ctl_name, ctl_files) in codeml.iteritems():
        # This dictionary also contains the "ExcludeSamples" option
        if ctl_name in ("ExcludeSamples", "SubsetRegions"):
            continue

        for (sequence, filename) in fasta_files.iteritems():
            output_tar = os.path.join(destination, "%s.%s.tar.gz" % (sequence, ctl_name))
            ctl_file = ctl_files["ControlFile"].format(Name=sequence)
            tree_file = ctl_files["TreeFile"].format(Name=sequence)

            node = CodemlNode(control_file=ctl_file,
                              trees_file=tree_file,
                              sequence_file=filename,
                              output_tar=output_tar,
                              exclude_groups=codeml["ExcludeSamples"],
                              dependencies=dependencies)
            codeml_nodes.append(node)

    return codeml_nodes


def chain_codeml(_pipeline, options, makefiles):
    destination = options.destination  # Move to makefile
    for makefile in makefiles:
        nodes = []
        filtering = makefile["Project"]["FilterSingletons"]
        options.destination = os.path.join(destination, makefile["Project"]["Title"])

        for regions in makefile["Project"]["Regions"].itervalues():
            nodes.extend(build_codeml_nodes(options, makefile, regions, filtering, makefile["Nodes"]))
        makefile["Nodes"] = tuple(nodes)
    options.destination = destination
