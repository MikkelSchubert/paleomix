import os

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd, ParallelCmds


class GenotypeNode(CommandNode):
    pileup_args = "-EA"
    caller_args = "-g"

    def __init__(self, config, reference, infile, outfile, regions = None, dependencies = ()):
        assert outfile.lower().endswith(".vcf.bgz")

        call_pileup = ["samtools", "mpileup", 
                       GenotypeNode.pileup_args,
                       "-uf", "%(IN_REFERENCE)s",
                       "%(IN_BAMFILE)s"]
        if regions:
            call_pileup[-1:-1] = ["-l", "%(IN_REGIONS)s"]

        pileup   = AtomicCmd(call_pileup,
                             IN_REFERENCE = reference,
                             IN_BAMFILE   = infile,
                             IN_REGIONS   = regions,
                             stdout       = AtomicCmd.PIPE,
                             stderr       = outfile + ".pileup_log")
        
        genotype = AtomicCmd(["bcftools", "view", GenotypeNode.caller_args, "-"],
                             stdin        = pileup,
                             stdout       = AtomicCmd.PIPE,
                             stderr       = outfile + ".genotype_log")

        bgzip    = AtomicCmd(["bgzip"],
                             stdin        = genotype,
                             stdout       = outfile)

        description = "<Genotyper: '%s' -> '%s'>" % (infile, outfile)
        CommandNode.__init__(self, 
                             description  = description,
                             command      = ParallelCmds([pileup, genotype, bgzip]),
                             dependencies = dependencies)



class TabixIndexNode(CommandNode):
    def __init__(self, config, infile, preset = "vcf", dependencies = ()):
        assert infile.lower().endswith(".vcf.bgz")

        self._infile = infile
        cmd_tabix = AtomicCmd(["tabix", "-p", preset, "%(TEMP_IN_VCFFILE)s"],
                              TEMP_IN_VCFFILE = os.path.basename(infile),
                              IN_VCFFILE      = infile,
                              OUT_TBI         = infile + ".tbi")

        CommandNode.__init__(self, 
                             description  = "<TabixIndex (%s): '%s'>" % (preset, infile,),
                             command      = cmd_tabix,
                             dependencies = dependencies)

    def _setup(self, config, temp):
        infile  = os.path.abspath(self._infile)
        outfile = os.path.join(temp, os.path.basename(self._infile))
        os.symlink(infile, outfile)

        return CommandNode._setup(self, config, temp)


    def _teardown(self, config, temp):
        os.remove(os.path.join(temp, os.path.basename(self._infile)))

        return CommandNode._teardown(self, config, temp)


class FastaIndexNode(CommandNode):
    def __init__(self, config, infile, dependencies = ()):
        self._infile = infile
        cmd_faidx = AtomicCmd(["samtools", "faidx", "%(TEMP_IN_FASTA)s"],
                              TEMP_IN_FASTA = os.path.basename(infile),
                              IN_FASTA      = infile,
                              OUT_TBI       = infile + ".fai")

        CommandNode.__init__(self, 
                             description  = "<FastaIndex: '%s'>" % (infile,),
                             command      = cmd_faidx,
                             dependencies = dependencies)

    def _setup(self, _config, temp):
        infile  = os.path.abspath(self._infile)
        outfile = os.path.join(temp, os.path.basename(self._infile))
        os.symlink(infile, outfile)

        return CommandNode._setup(self, config, temp)


    def _teardown(self, _config, temp):
        os.remove(os.path.join(temp, os.path.basename(self._infile)))

        return CommandNode._teardown(self, config, temp)
