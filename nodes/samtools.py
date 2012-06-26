import os

import node
import fileutils

from atomiccmd import AtomicCmd, ParallelCmds


class GenotypeNode(node.Node):
    pileup_args = "-EA"
    caller_args = "-g"

    def __init__(self, config, destination, reference, infile, outfile, regions = None, dependencies = ()):
        assert outfile.lower().endswith(".vcf.bgz")

        call_pileup = ["samtools", "mpileup", 
                       GenotypeNode.pileup_args,
                       "-uf", "%(IN_REFERENCE)s",
                       "%(IN_BAMFILE)s"]
        if regions:
            call_pileup[-1:-1] = ["-l", "%(IN_REGIONS)s"]

        pileup   = AtomicCmd(destination,
                             call_pileup,
                             IN_REFERENCE = reference,
                             IN_BAMFILE   = infile,
                             IN_REGIONS   = regions,
                             stdout       = AtomicCmd.PIPE,
                             stderr       = outfile + ".pileup_log")
        
        genotype = AtomicCmd(destination,
                             ["bcftools", "view", GenotypeNode.caller_args, "-"],
                             stdin        = pileup,
                             stdout       = AtomicCmd.PIPE,
                             stderr       = outfile + ".genotype_log")

        bgzip    = AtomicCmd(destination,
                             ["bgzip"],
                             stdin        = genotype,
                             stdout       = outfile)

        description = "<Genotyper: '%s' -> '%s'>" \
            % (infile, os.path.join(destination, outfile))

        node.Node.__init__(self, 
                           description  = description,
                           command      = ParallelCmds([pileup, genotype, bgzip]),
                           dependencies = dependencies)



# FIXME: Should use temp folder ...
class TabixIndexNode(node.Node):
    def __init__(self, config, destination, infile, preset = "vcf", dependencies = ()):
        assert infile.lower().endswith(".vcf.bgz")

        self._infile = infile
        cmd_tabix = AtomicCmd(destination,
                              ["tabix", "-p", preset, "%(IN_VCFFILE)s"],
                              IN_VCFFILE = infile)

        node.Node.__init__(self, 
                           description  = "<TabixIndex (%s): '%s'>" % (preset, infile,),
                           command      = cmd_tabix,
                           dependencies = dependencies)

    def output_exists(self):
        filename = self._infile + ".tbi"
        
        return fileutils.valid_file(filename)
