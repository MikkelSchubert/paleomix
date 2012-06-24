import os

import pypeline
import fileutils

import nodes.picard
import nodes.gatk
import nodes.samtools
import nodes.bam_statistics


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

        validate = nodes.picard.ValidateBAMNode(config,
                                                bamfile     = infile, 
                                                ignore      = ["MATE_NOT_FOUND"])

        stats    = nodes.bam_statistics.PairedStatistics(config,
                                                         infile      = infile)
    
        realigned_bam = fileutils.add_postfix(os.path.split(infile)[-1], ".realigned")
        realigner = nodes.gatk.add_indel_realigner(config,
                                                   destination  = "dest",
                                                   reference    = "data/reference.fasta",
                                                   infile       = infile,
                                                   outfile      = realigned_bam,
                                                   dependencies = [validate, stats])

        validate = nodes.picard.ValidateBAMNode(config, 
                                                bamfile      = os.path.join("dest", realigned_bam),
                                                ignore       = ["MATE_NOT_FOUND"],
                                                dependencies = [realigner])
    
        vcffile = fileutils.swap_ext(realigned_bam, ".vcf.bgz")
        genotype = nodes.samtools.Genotype(config       = config,
                                           destination  = "dest",
                                           reference    = "data/reference.fasta",
                                           regions      = "data/intervals.bed",
                                           infile       = os.path.join("dest", realigned_bam),
                                           outfile      = vcffile,
                                           dependencies = [validate])

        tabix = nodes.samtools.TabixIndex(config       = config,
                                          destination  = "dest",
                                          infile       = os.path.join("dest", vcffile),
                                          dependencies = [genotype])

        pipeline.add_node(tabix)

    pipeline.run()
