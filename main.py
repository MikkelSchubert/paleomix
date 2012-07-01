import os

import pypeline
import fileutils

from nodes.picard import ValidateBAMNode
from nodes.gatk import IndelRealignerNode
from nodes.bedtools import SlopBedNode
from nodes.samtools import GenotypeNode, TabixIndexNode
from nodes.bam_statistics import PairedStatisticsNode


if __name__ == '__main__':
    class Config:
        jar_root  = "opt/picard-tools"
        gatk_root = "opt/GATK"
        temp_root = "temp"


    pipeline = pypeline.Pypeline(Config)
    for infile in ('data/small.bam', 'data/medium.bam', 'data/merged.bam'):
        filename = fileutils.add_postfix(os.path.basename(infile), ".realigned")
        realigned_bam = os.path.join("dest", filename)

        validate   = ValidateBAMNode(Config,
                                     bamfile     = infile, 
                                     ignore      = ["MATE_NOT_FOUND"])

        statistics = PairedStatisticsNode(Config,
                                          infile      = infile)
    
        realigner = IndelRealignerNode(Config,
                                       reference    = "data/reference.fasta",
                                       infile       = infile,
                                       outfile      = realigned_bam,
                                       dependencies = [validate, statistics])

        validate = ValidateBAMNode(Config,
                                   bamfile      = realigned_bam,
                                   ignore       = ["MATE_NOT_FOUND"],
                                   dependencies = [realigner])


        padding = 10
        vcfraw  = "data/intervals.bed"
        vcfslop = fileutils.swap_ext(realigned_bam, ".vcf.bgz.bed_slop")
        vcffile = fileutils.swap_ext(realigned_bam, ".vcf.bgz")
        
        slop_bed = SlopBedNode(config        = Config,
                               genome        = "data/reference.genome",
                               infile        = vcfraw,
                               outfile       = vcfslop,
                               from_start    = 10,
                               from_end      = 10)

        genotype = GenotypeNode(config       = Config,
                                reference    = "data/reference.fasta",
                                regions      = vcfslop,
                                infile       = realigned_bam,
                                outfile      = vcffile,
                                dependencies = [validate, slop_bed])

        tabix = TabixIndexNode(config       = Config,
                               infile       = vcffile,
                               dependencies = [genotype])

        pipeline.add_node(tabix)

    pipeline.run()
