#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 Mikkel Schubert <MSchubert@snm.ku.dk>
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
"""Wrapper script around PALEOMIX commands.

This script takes care of checking that various requirements are met, that the
PALEOMIX module ('paleomix') is available, and forwards arguments to the
appropriate commands.
"""
import os
import sys
import textwrap


# List of tuples of commands: (name, module, help string).
# If module is None, the command cannot be invoked directly (e.g. help), and
# if help string is none, the command is considered a help header.
def _commands():
    yield ("Pipelines", None, None)
    yield ("bam_pipeline", "paleomix.tools.bam_pipeline.pipeline",
           "Pipeline for trimming and mapping of NGS reads.")
    yield ("trim_pipeline", "paleomix.tools.bam_pipeline.trim_pipeline",
           "Equivalent to 'bam_pipeline', but only runs the trimming steps.")
    yield ("phylo_pipeline", "paleomix.tools.phylo_pipeline.pipeline",
           "Pipeline for genotyping and phylogenetic inference from BAMs.")
    yield ("zonkey", "paleomix.tools.zonkey.pipeline",
           "Pipeline for detecting F1 (equine) hybrids.")

    # Currently not documented; used internally by Zonkey
    yield ("zonkey_db", "paleomix.tools.zonkey.build_db", None)
    yield ("zonkey_tped", "paleomix.tools.zonkey.build_tped", None)
    yield ("zonkey_mito", "paleomix.tools.zonkey.build_mito", None)

    # In development: Integration of metabit into PALEOMIX
    yield ("metabit", "paleomix.tools.metabit.metabit", None)

    yield ("BAM/SAM tools", None, None)
    yield ("cleanup", "paleomix.tools.cleanup",
           "Reads SAM file from STDIN, and outputs sorted, tagged, and filter "
           "BAM, for which NM and MD tags have been updated.")
    yield ("coverage", "paleomix.tools.coverage",
           "Calculate coverage across reference sequences or regions of "
           "interest.")
    yield ("depths", "paleomix.tools.depths",
           "Calculate depth histograms across reference sequences or regions "
           "of interest.")
    yield ("duphist", "paleomix.tools.duphist",
           "Generates PCR duplicate histogram; used with the 'Preseq' tool.")
    yield ("rmdup_collapsed", "paleomix.tools.rmdup_collapsed",
           "Filters PCR duplicates for collapsed paired-ended reads generated "
           "by the AdapterRemoval tool.")

    yield ("VCF/GTF/BED/Pileup tools", None, None)
    yield ("genotype", "paleomix.tools.genotype",
           "Creates bgzipped VCF for a set of (sparse) BED regions, or for "
           "entire chromosomes / contigs using SAMTools / BCFTools.")
    yield ("gtf_to_bed", "paleomix.tools.gtf_to_bed",
           "Convert GTF file to BED files grouped by feature "
           "(coding, RNA, etc).")
    yield ("sample_pileup", "paleomix.tools.sample_pileup",
           "Randomly sample sites in a pileup to generate a FASTA sequence.")
    yield ("vcf_filter", "paleomix.tools.vcf_filter",
           "Quality filters for VCF records, similar to "
           "'vcfutils.pl varFilter'.")
    yield ("vcf_to_fasta", "paleomix.tools.vcf_to_fasta",
           "Create most likely FASTA sequence from tabix-indexed VCF file.")

    yield ("Misc tools", None, None)
    yield ("cat", "paleomix.tools.cat",
           "Generalized cat command for gz, bz2 and uncompressed files.")

    # In development:
    #   Prepares FASTQ reads recorded in BAM pipeline makefiles
    #   for submission to the European Nucleotide Archive.
    yield ("ena", "paleomix.tools.ena", None)

# Error message shown if the Pysam module ('pysam') cannot be imported
_IMPORT_ERROR_PYSAM = """
Error importing required python module 'pysam':
    - %s

The module may be installed for the current user using 'pip':
    $ pip install --user pysam

Alternatively, download the latest version from the Pysam repository at GitHub:
    - https://github.com/pysam-developers/pysam

A local install may be performed using the following command:
     $ python setup.py install --user
"""

# Error message in case it is not possible to import the PALEOMIX module itself
_IMPORT_ERROR_PALEOMIX = """
Error importing PALEOMIX module 'paleomix':
    - %s

Please make sure that PALEOMIX is correctly installed, and that the PYTHONPATH
environmental variable points to the location of the 'paleomix' module.
"""

_INCONSISTENT_IMPORT_ERROR = """
Inconsistency importing PALEOMIX module 'paleomix'; the currently running
script is not located in the same folder as the 'paleomix' module. This
suggests that you have multiple, conflicting copies of PALEOMIX installed!

  - The running script:    %r
  - The 'paleomix' module: %r

It is strongly suggested that you remove all installed copies of PALEOMIX,
and perform a clean install. If this is not possible, the 'virtualenv' tool
for Python may be used to prevent conflict between the installed versions.
"""


_PALEOMIX_CITATION = """

If you make use of PALEOMIX in your work, please cite
  Schubert et al, "Characterization of ancient and modern genomes by SNP
  detection and phylogenomic and metagenomic analysis using PALEOMIX".
  Nature Protocols. 2014 May; 9(5): 1056-82. doi: 10.1038/nprot.2014.063
"""


def _are_requirements_met():
    """Checks the current Python version, that the PALEOMIX modules are
    available, and modules required by the pipeline (Pysam) is available and
    up-to-date.
    """
    if tuple(sys.version_info)[:2] != (2, 7):
        sys.stderr.write("ERROR: PALEOMIX requires Python version 2.7.x.\n")
        sys.stderr.write("However, the current version of python is\n\tv%s\n\n"
                         % (sys.version.replace("\n", "\n\t"),))
        sys.stderr.write("Please install Python v2.7 to continue.\n")
        return False

    modules = [('pysam', _IMPORT_ERROR_PYSAM),
               ('paleomix', _IMPORT_ERROR_PALEOMIX)]

    for (module, message) in modules:
        try:
            __import__(module)
        except ImportError:
            error = sys.exc_info()[1]  # Python 2/3 compatible exception syntax
            sys.stderr.write(message % (error,))
            return False

    # Sanity check, to catch multiple, conflicting PALEOMIX installations
    import paleomix
    if not os.path.samefile(os.path.dirname(__file__),
                            os.path.dirname(paleomix.__file__)):
        sys.stderr.write(_INCONSISTENT_IMPORT_ERROR
                         % (os.path.dirname(__file__),
                            os.path.dirname(paleomix.__file__)))
        return False

    import pysam
    version = [int(field) for field in pysam.__version__.split(".")]
    if version[:3] < [0, 8, 3]:
        error = "Pysam is outdated (v%s), version must be at least v0.8.3!"
        error %= (pysam.__version__,)
        sys.stderr.write(_IMPORT_ERROR_PYSAM % (error,))
        return False

    return True


def _print_help():
    """Prints description of commands and reference to PALEOMIX paper."""
    import paleomix

    template = "    paleomix %s%s-- %s\n"
    max_len = max(len(key) for (key, module, _) in _commands() if module)
    help_len = 80 - len(template % (" " * max_len, " ", ""))
    help_padding = (80 - help_len) * " "

    sys.stderr.write("PALEOMIX - pipelines and tools for NGS data analyses.\n")
    sys.stderr.write("Version: %s\n\n" % (paleomix.__version__,))
    sys.stderr.write("Usage: paleomix <command> [options]\n")
    for (key, module, help_str) in _commands():
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
    """Main function; takes a list of arguments excluding argv[0]."""
    if not _are_requirements_met():
        return 1

    # Process name defaults to the name of the python executable
    import paleomix.common.system
    paleomix.common.system.set_procname("paleomix")

    if not argv or argv[0] == "help":
        _print_help()
        return 0

    command = argv[0]
    for (cmd_name, cmd_module, _) in _commands():
        if cmd_module and (command == cmd_name):
            module = __import__(cmd_module, globals(), locals(), ["main"], 0)
            return module.main(argv[1:])

    sys.stderr.write("ERROR: Unknown PALEOMIX command %r!\n" % (command,))
    return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
