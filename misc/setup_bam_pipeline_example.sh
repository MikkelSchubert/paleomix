#!/bin/bash
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

SP_SEED=${RANDOM}

rm -rv 000_data
mkdir -p 000_data
for barcode in ACGATA GCTCTG TGCTCA;
do
    python $(dirname $0)/synthesize_reads.py 000_prefixes/rCRS.fasta 000_data/ \
	--library-barcode=${barcode} \
	--specimen-seed=${SP_SEED} \
	--lanes-reads-mu=2500 \
	--lanes-reads-sigma=500 \
	--lanes-reads-per-file=1000 \
	--lanes=2 \
	--damage
done

rm -v 000_data/GCTCTG_L*R2*.gz
rm -v 000_data/TGCTCA_L1_R2*.gz

bam_pipeline run $(dirname $0)/setup_bam_pipeline_example.makefile.yaml --destination .

mkdir -p 000_data/ACGATA_L2/
mv ExampleProject/reads/Synthetic_Sample_1/ACGATA/Lane_2/reads.singleton.truncated.gz 000_data/ACGATA_L2/reads.singleton.truncated.gz
mv ExampleProject/reads/Synthetic_Sample_1/ACGATA/Lane_2/reads.collapsed.gz 000_data/ACGATA_L2/reads.collapsed.gz
mv ExampleProject/reads/Synthetic_Sample_1/ACGATA/Lane_2/reads.collapsed.truncated.gz 000_data/ACGATA_L2/reads.collapsed.truncated.gz

mv ExampleProject/rCRS/Synthetic_Sample_1/GCTCTG/Lane_2/single.minQ0.bam 000_data/GCTCTG_L2.bam

rm -v 000_data/ACGATA_L2_R*.fastq.gz
rm -v 000_data/GCTCTG_L2_R1_*.fastq.gz
rm -rv ExampleProject
