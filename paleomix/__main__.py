#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
#
from __future__ import annotations

import faulthandler
import logging
import signal
import sys

import pysam
import setproctitle

import paleomix
import paleomix.common.logging

COMMANDS = {
    # CMBR NGS pipeline
    "ngs": "paleomix.pipelines.ngs",
    "ngs:finalize_bam": "paleomix.pipelines.ngs.tools.finalize_bam",
    # BAM/FASTQ pipeline
    "bam": "paleomix.pipelines.bam.main",
    "bam_pipeline": "paleomix.pipelines.bam.main",
    "trim": "paleomix.pipelines.bam.trim_pipeline",
    "trim_pipeline": "paleomix.pipelines.bam.trim_pipeline",
    "worker": "paleomix.tools.worker",
    # Zonkey
    "zonkey": "paleomix.pipelines.zonkey.pipeline",
    "zonkey:db": "paleomix.pipelines.zonkey.build_db",
    "zonkey:mito": "paleomix.pipelines.zonkey.build_mito",
    "zonkey:tped": "paleomix.pipelines.zonkey.build_tped",
    # BAM file tools
    "cleanup": "paleomix.tools.cleanup",
    "coverage": "paleomix.tools.coverage",
    "depths": "paleomix.tools.depths",
    "dupcheck": "paleomix.tools.dupcheck",
    # VCF/etc. tools
    "rmdup_collapsed": "paleomix.tools.rmdup_collapsed",
    # Misc tools
    ":rscript": "paleomix.tools.rscript",
    ":validate_fastq": "paleomix.tools.validate_fastq",
    ":validate_fasta": "paleomix.tools.validate_fasta",
}


_HELP = """PALEOMIX - pipelines and tools for NGS data analyses
Version: {version}

Pipelines:
    paleomix bam              -- Pipeline for trimming and mapping of NGS reads.
    paleomix trim             -- Equivalent to the 'bam' pipeline, but only runs
                                 the FASTQ trimming steps.
    paleomix zonkey           -- Pipeline for detecting F1 (equine) hybrids.
    paleomix worker           -- Distribute pipelines across multiple systems.

BAM/SAM tools:
    paleomix coverage         -- Calculate coverage across reference sequences
                                 or regions of interest.
    paleomix depths           -- Calculate depth histograms across reference
                                 sequences or regions of interest.
    paleomix rmdup_collapsed  -- Filters PCR duplicates for collapsed paired-
                                 ended reads generated by the AdapterRemoval
                                 tool.

If you make use of PALEOMIX in your work, please cite
  Schubert et al, "Characterization of ancient and modern genomes by SNP
  detection and phylogenomic and metagenomic analysis using PALEOMIX".
  Nature Protocols. 2014 May; 9(5): 1056-82. doi: 10.1038/nprot.2014.063
"""


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    if not argv or argv[0] in ("-h", "--help", "help"):
        print(_HELP.format(version=paleomix.__version__))
        return 0
    elif argv[0] in ("--version",):
        print(f"paleomix v{paleomix.__version__}")
        return 0

    # Change process name to 'paleomix' regardless of how the script was invoked
    setproctitle.setproctitle("paleomix " + " ".join(argv[:1]))
    # Setup basic logging to STDERR
    paleomix.common.logging.initialize_console_logging()
    # Silence log-messages from HTSLIB
    pysam.set_verbosity(0)  # pyright: ignore[reportUnknownMemberType]
    # Add handlers to print backtraces on SIGSEGV, SIGFPE, SIGABRT, SIGBUS and SIGILL
    faulthandler.enable()
    # Add handler to print backtraces on SIGUSR1
    faulthandler.register(signal.SIGUSR1)

    command = COMMANDS.get(argv[0])
    if command is None:
        log = logging.getLogger(__name__)
        log.error("Unknown command %r", argv[0])
        return 1

    module = __import__(command, fromlist=["main"])

    return module.main(argv[1:])


if __name__ == "__main__":
    sys.exit(main())
