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

rm -rvf alignment/000_reads
mkdir -p alignment/000_reads
for PREFIX in `ls alignment/000_prefixes/*.fasta | grep -v rCRS`;
do
    SP_SEED=${RANDOM}
    NAME=$(echo ${PREFIX} | sed -e's#alignment/000_prefixes/##' -e's#\..*##')
    mkdir -p alignment/000_reads/${NAME/*\//}/

    ./synthesize_reads.py ${PREFIX} alignment/000_reads/${NAME}/ \
	--specimen-seed=${SP_SEED} \
	--lanes-reads-mu=50000 \
	--lanes-reads-sigma=500 \
	--lanes-reads-per-file=10000 \
	--reads-len=50 \
	--lanes=1
done

# These links would not survive the package installation, so setup here
ln -sf ../../alignment/000_prefixes/ phylogeny/data/prefixes
ln -sf ../../alignment phylogeny/data/samples

# Create link to reference sequence
mkdir -p phylogeny/data/refseqs
ln -sf ../../../alignment/000_prefixes/rCRS.fasta phylogeny/data/refseqs/rCRS.rCRS.fasta
