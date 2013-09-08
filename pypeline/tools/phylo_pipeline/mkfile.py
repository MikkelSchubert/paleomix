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
import sys


_TEMPLATE = """# -*- mode: Yaml; -*-
Project:
  Title: PROJECT_NAME

  Taxa:
    <GROUP>:
      <SUBGROUP>:
        TAXA_NAME:
          Species Name: ...
          Common Name:  ...
          Gender:       ...


  Intervals:
     INTERVALS_NAME:
       Genome: GENOME_NAME
       Protein coding: no
       # Do not include indels in final sequence
       # Indels are still called, and used to filter SNPs
       Include indels: yes
       Homozygous Contigs:
         GENDER_X: []


#  Filter sites in a taxa, excluding any nucleotide not
#  observed in the specified list of taxa or groups.
#  Filter Singletons: {}
#    <TAXA> : [<TAXA_1>, ...]


Genotyping:
  # Padding used for genotyping, to ensure that we call adjacent indels
  Padding: 10

  Random:
    # Min distance of variants to indels
    --min-distance-to-indels: 2

  MPileup:
    -E: # extended BAQ for higher sensitivity but lower specificity
    -A: # count anomalous read pairs

  BCFTools:
    -g: # Call genotypes at variant sites

  VCF_Filter:
    # Maximum coverage acceptable for genotyping calls
    # If zero, the default vcf_filter value is used
    MaxReadDepth: 0

    # Minimum coverage acceptable for genotyping calls
    --min-read-depth: 8
    # Min RMS mapping quality
    --min-mapping-quality: 10
    # Min QUAL score (Phred) for genotyping calls
    --min-quality: 30
    # Min distance of variants to indels
    --min-distance-to-indels: 2
    # Min distance between indels
    --min-distance-between-indels: 10
    # Min P-value for strand bias (given PV4)
    --min-strand-bias: 1e-4
    # Min P-value for baseQ bias (given PV4)
    --min-baseq-bias: 1e-4
    # Min P-value for mapQ bias (given PV4)
    --min-mapq-bias: 1e-4
    # Min P-value for end distance bias (given PV4)
    --min-end-distance-bias: 1e-4
    # Max frequency of the major allele at heterozygous sites
    --min-allele-frequency: 0.2


MSAlignment:
  Enabled: yes

  MAFFT:
    Algorithm: G-INS-i



Phylogenetic Inference:
  ExcludeGroups: []

  ExaML:
    Threads: 8
    # Number of times to perform full phylogenetic inference
    Replicates: 1
    # Number of bootstraps to compute
    Bootstraps: 100
    Model: GAMMA


PAML:
  codeml:
    ExcludeGroups: []

    # Allow auto-generation of path from options.destination, Project/Title
    Control Files:
       "Null": "makefiles/Common.codeml.null.ctl"
       "Test": "makefiles/Common.codeml.test.ctl"
    Tree File: "makefiles/3Cabs.codeml.trees"
"""


def main(_argv):
    print _TEMPLATE

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
