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
import sys
import datetime

import pypeline.ui as ui
from pypeline.common.text import padded_table


def main(argv):
    print """# -*- mode: Yaml; -*-
Project:
  Title: Horses

  Taxa:
    <GROUP>:
      <SUBGROUP>:
        SPECIES NAME:
          Species Name: ...
          Common Name:  ...
          Gender:       ...


  Intervals:
     Ensembl.v71.EquCab.500_known.protein_coding.CDS:
       Genome: Equus_cab_nucl_wChrUn
       Protein coding: yes
       Homozygous Contigs:
         Female: []
         Male: [chrX]

  Filter Singletons: {}


Genotyping:
  # Default genotyping method
  Default: SAMTools

  Random:
    # Padding used for genotyping, to ensure that we call adjacent indels
    --padding: 10
    # Min distance of variants to indels
    --min-distance-to-indels: 2

  MPileup:
    -E: # extended BAQ for higher sensitivity but lower specificity
    -A: # count anomalous read pairs

  BCFTools:
    -g: # Call genotypes at variant sites

  VCF_Filter:
    # Padding used for genotyping, to ensure that we call adjacent indels
    --padding:   10
    # Minimum coverage acceptable for genotyping calls
    --min-read-depth: 8
    # Minimum coverage acceptable for genotyping calls
    MaxDepth: 100
    # Min QUAL score (Phred) for genotyping calls
    --min-quality: 30
    # Min distance of variants to indels
    --min-distance-to-indels: 2



MSAlignment:
  Default: MAFFT

  MAFFT:
    Algorithm: G-INS-i



Phylogenetic Inference:
  Default: ExaML

  ExaML:
    # Number of times to perform full phylogenetic inference
    Replicates: 1
    # Number of bootstraps to compute
    Bootstraps: 100
    Model: GAMMA


PAML:
  codeml:
    # Allow auto-generation of path from options.destination, Project/Title
    Control Files:
      - "results/Horses.codeml.ctl"
    Tree File:    "rusults/Horses.codeml.trees"

"""

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
