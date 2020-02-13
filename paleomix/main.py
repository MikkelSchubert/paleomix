#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import logging
import sys
import textwrap

import paleomix.common.system
import paleomix.logger


# List of tuples of commands: (name, module, help string). If module is None,
# then the ntry is treated as a header, and if the help string is empty then
# the command is not listed when invoking "help".
_PALEOMIX_COMMANDS = (
    ("Pipelines", None, None),
    (
        "bam_pipeline",
        "paleomix.pipelines.ngs.pipeline",
        "Pipeline for trimming and mapping of NGS reads.",
    ),
    (
        "trim_pipeline",
        "paleomix.pipelines.ngs.trim_pipeline",
        "Equivalent to 'bam_pipeline', but only runs the trimming steps.",
    ),
    (
        "phylo_pipeline",
        "paleomix.pipelines.phylo.pipeline",
        "Pipeline for genotyping and phylogenetic inference from BAMs.",
    ),
    (
        "zonkey",
        "paleomix.pipelines.zonkey.pipeline",
        "Pipeline for detecting F1 (equine) hybrids.",
    ),
    # Currently not documented; used internally by Zonkey
    ("zonkey_db", "paleomix.pipelines.zonkey.build_db", None),
    ("zonkey_tped", "paleomix.pipelines.zonkey.build_tped", None),
    ("zonkey_mito", "paleomix.pipelines.zonkey.build_mito", None),
    ("BAM/SAM tools", None, None),
    (
        "dupcheck",
        "paleomix.tools.dupcheck",
        "Identifies potential duplicate data in sorted BAM files, defined"
        "as reads aligned to the same position, with the same name, "
        "sequence, and qualities.",
    ),
    (
        "cleanup",
        "paleomix.tools.cleanup",
        "Reads SAM file from STDIN, and outputs sorted, tagged, and filter "
        "BAM, for which NM and MD tags have been updated.",
    ),
    (
        "coverage",
        "paleomix.tools.coverage",
        "Calculate coverage across reference sequences or regions of " "interest.",
    ),
    (
        "depths",
        "paleomix.tools.depths",
        "Calculate depth histograms across reference sequences or regions "
        "of interest.",
    ),
    (
        "duphist",
        "paleomix.tools.duphist",
        "Generates PCR duplicate histogram; used with the 'Preseq' tool.",
    ),
    (
        "rmdup_collapsed",
        "paleomix.tools.rmdup_collapsed",
        "Filters PCR duplicates for collapsed paired-ended reads generated "
        "by the AdapterRemoval tool.",
    ),
    ("VCF/GTF/BED/Pileup tools", None, None),
    (
        "gtf_to_bed",
        "paleomix.tools.gtf_to_bed",
        "Convert GTF file to BED files grouped by feature " "(coding, RNA, etc).",
    ),
    (
        "vcf_filter",
        "paleomix.tools.vcf_filter",
        "Quality filters for VCF records, similar to " "'vcfutils.pl varFilter'.",
    ),
    (
        "vcf_to_fasta",
        "paleomix.tools.vcf_to_fasta",
        "Create most likely FASTA sequence from tabix-indexed VCF file.",
    ),
    ("Misc tools", None, None),
    (
        "cat",
        "paleomix.tools.cat",
        "Generalized cat command for gz, bz2 and uncompressed files.",
    ),
    (
        "retable",
        "paleomix.tools.retable",
        "Pretty print whitespace separated tabular data.",
    ),
)


_PALEOMIX_CITATION = """

If you make use of PALEOMIX in your work, please cite
  Schubert et al, "Characterization of ancient and modern genomes by SNP
  detection and phylogenomic and metagenomic analysis using PALEOMIX".
  Nature Protocols. 2014 May; 9(5): 1056-82. doi: 10.1038/nprot.2014.063
"""


def _print_help():
    """Prints description of commands and reference to PALEOMIX paper."""
    import paleomix

    template = "    paleomix %s%s-- %s\n"
    max_len = max(len(key) for (key, module, _) in _PALEOMIX_COMMANDS if module)
    help_len = 80 - len(template % (" " * max_len, " ", ""))
    help_padding = (80 - help_len) * " "

    sys.stderr.write("PALEOMIX - pipelines and tools for NGS data analyses.\n")
    sys.stderr.write("Version: %s\n\n" % (paleomix.__version__,))
    sys.stderr.write("Usage: paleomix <command> [options]\n")
    for (key, module, help_str) in _PALEOMIX_COMMANDS:
        if help_str is None:
            if module is None:
                sys.stderr.write("\n%s:\n" % (key,))
        else:
            lines = textwrap.wrap(help_str, help_len)
            padding = (max_len - len(key) + 2) * " "
            sys.stderr.write(template % (key, padding, lines[0]))

            for line in lines[1:]:
                sys.stderr.write("%s%s\n" % (help_padding, line))

    sys.stderr.write(_PALEOMIX_CITATION)


def main(argv):
    # Change process name from 'python' to 'paleomix'
    paleomix.common.system.set_procname("paleomix")
    # Setup basic logging to STDERR
    paleomix.logger.initialize_console_logging()

    if not argv or argv[0] == "help":
        _print_help()
        return 0

    command = argv[0]
    for (cmd_name, cmd_module, _) in _PALEOMIX_COMMANDS:
        if cmd_module and (command == cmd_name):
            module = __import__(cmd_module, globals(), locals(), ["main"], 0)
            return module.main(argv[1:])

    log = logging.getLogger(__name__)
    log.error("Unknown command %r", command)
    return 1


def entry_point():
    return main(sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
