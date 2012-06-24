import os

import pypeline
import fileutils

from nodes.picard import ValidateBAMNode
from nodes.gatk import IndelRealignerNode
from nodes.samtools import GenotypeNode, TabixIndexNode
from nodes.bam_statistics import PairedStatisticsNode


if __name__ == '__main__':
    fileutils.TEMP_ROOT = "temp/"

    class config:
        def __init__(self):
            pass
    config.jar_root = "opt/picard-tools"
    config.gatk_root = "opt/GATK"


    pipeline = pypeline.Pypeline(config)
    for infile in ('data/small.bam', 'data/medium.bam'):
        filename = os.path.split(infile)[-1]

        validate   = ValidateBAMNode(config,
                                     bamfile     = infile, 
                                     ignore      = ["MATE_NOT_FOUND"])

        statistics = PairedStatisticsNode(config,
                                          infile      = infile)
    
        realigned_bam = fileutils.add_postfix(os.path.basename(infile), ".realigned")
        realigner = IndelRealignerNode(config,
                                       destination  = "dest",
                                       reference    = "data/reference.fasta",
                                       infile       = infile,
                                       outfile      = realigned_bam,
                                       dependencies = [validate, statistics])

        validate = ValidateBAMNode(config, 
                                   bamfile      = os.path.join("dest", realigned_bam),
                                   ignore       = ["MATE_NOT_FOUND"],
                                   dependencies = [realigner])
        
        vcffile = fileutils.swap_ext(realigned_bam, ".vcf.bgz")
        genotype = GenotypeNode(config       = config,
                                destination  = "dest",
                                reference    = "data/reference.fasta",
                                regions      = "data/intervals.bed",
                                infile       = os.path.join("dest", realigned_bam),
                                outfile      = vcffile,
                                dependencies = [validate])

        tabix = TabixIndexNode(config       = config,
                               destination  = "dest",
                               infile       = os.path.join("dest", vcffile),
                               dependencies = [genotype])

        pipeline.add_node(tabix)

    pipeline.run()
