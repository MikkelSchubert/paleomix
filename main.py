import os

import pypeline
import fileutils

from nodes.picard import ValidateBAMNode
from nodes.gatk import IndelRealignerNode
from nodes.samtools import GenotypeNode, TabixIndexNode
from nodes.bam_statistics import PairedStatisticsNode


if __name__ == '__main__':
    class Config:
	jar_root  = "opt/picard-tools"
    	gatk_root = "opt/GATK"
	temp_root = "temp"


    pipeline = pypeline.Pypeline(Config)
    for infile in ('data/small.bam', 'data/medium.bam'):
        filename = os.path.split(infile)[-1]

        validate   = ValidateBAMNode(Config,
                                     bamfile     = infile, 
                                     ignore      = ["MATE_NOT_FOUND"])

        statistics = PairedStatisticsNode(Config,
                                          infile      = infile)
    
        realigned_bam = fileutils.add_postfix(os.path.basename(infile), ".realigned")
        realigner = IndelRealignerNode(Config,
                                       destination  = "dest",
                                       reference    = "data/reference.fasta",
                                       infile       = infile,
                                       outfile      = realigned_bam,
                                       dependencies = [validate, statistics])

        validate = ValidateBAMNode(Config, 
                                   bamfile      = os.path.join("dest", realigned_bam),
                                   ignore       = ["MATE_NOT_FOUND"],
                                   dependencies = [realigner])
        
        vcffile = fileutils.swap_ext(realigned_bam, ".vcf.bgz")
        genotype = GenotypeNode(config       = Config,
                                destination  = "dest",
                                reference    = "data/reference.fasta",
                                regions      = "data/intervals.bed",
                                infile       = os.path.join("dest", realigned_bam),
                                outfile      = vcffile,
                                dependencies = [validate])

        tabix = TabixIndexNode(config       = Config,
                               destination  = "dest",
                               infile       = os.path.join("dest", vcffile),
                               dependencies = [genotype])

        pipeline.add_node(tabix)

    pipeline.run()
