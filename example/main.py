#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
#
import os

import pypeline
import pypeline.fileutils as fileutils

from pypeline.nodes.picard import ValidateBAMNode
from pypeline.nodes.gatk import IndelRealignerNode
from pypeline.nodes.bedtools import SlopBedNode
from pypeline.nodes.samtools import GenotypeNode, TabixIndexNode
from pypeline.nodes.bam_statistics import PairedStatisticsNode


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
