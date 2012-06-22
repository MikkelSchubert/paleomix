import os

import pypeline
import fileutils

import nodes.picard
import nodes.gatk
import nodes.samtools

if __name__ == '__main__':
    fileutils.TEMP_ROOT = "temp/"

    class config:
        pass
    config.jar_root = "/home/mschubert/research/tools/bamPipeline/picard-tools-1.69"
    config.gatk_root = "/home/mschubert/research/opt/GATK/dist/"


    pipeline = pypeline.Pypeline(config)
    for infile in ('data/small.bam', 'data/medium.bam'):
        filename = os.path.split(infile)[-1]

        validate = nodes.picard.ValidateBAMNode(config,
                                                bamFile     = infile, 
                                                ignore      = ["MATE_NOT_FOUND"])
    
        realnBam = fileutils.add_postfix(os.path.split(infile)[-1], ".realigned")
        realigner = nodes.gatk.add_indel_realigner(config,
                                                   destination  = "dest",
                                                   reference    = "data/reference.fasta",
                                                   inFile       = infile,
                                                   outFile      = realnBam,
                                                   dependencies = [validate])

        validate = nodes.picard.ValidateBAMNode(config, 
                                                bamFile      = os.path.join("dest", realnBam),
                                                ignore       = ["MATE_NOT_FOUND"],
                                                dependencies = [realigner])
    
        genotype = nodes.samtools.Genotype(config       = config,
                                           destination  = "dest",
                                           reference    = "data/reference.fasta",
                                           infile       = os.path.join("dest", realnBam),
                                           outfile      = fileutils.swap_ext(realnBam, ".bcf"),
                                           dependencies = [validate])

        pipeline.add_node(genotype)

    pipeline.run()
