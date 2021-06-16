#!/bin/bash
#
# Copyright (c) 2013 Mikkel Schubert <MikkelSch@gmail.com>
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

# Exit on unset variables
set -o nounset

# Print debug message and terminate script on non-zero return codes
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# [1/2] Most versions of Bash read scripts line by line, causing misbehavior if
# the file changes during runtime. The {} forces Bash to read the entire thing
{
	rm -rvf alignment/reads
	mkdir -p alignment/reads
	for PREFIX in `ls alignment/prefixes/*.fasta | grep -v rCRS`;
	do
		SP_SEED=${RANDOM}
		NAME=$(echo ${PREFIX} | sed -e's#alignment/prefixes/##' -e's#\..*##')
		mkdir -p alignment/reads/${NAME/*\//}/

		./synthesize_reads.py ${PREFIX} alignment/reads/${NAME}/ \
			--specimen-seed=${SP_SEED} \
			--lanes-reads-mu=50000 \
			--lanes-reads-sigma=500 \
			--lanes-reads-per-file=10000 \
			--reads-len=50 \
			--lanes=1
	done

	# These links would not survive the package installation, so setup here
	ln -sf ../../alignment/prefixes/ phylogeny/data/prefixes
	ln -sf ../../alignment phylogeny/data/samples

	# Create link to reference sequence
	mkdir -p phylogeny/data/refseqs
	ln -sf ../../../alignment/prefixes/rCRS.fasta phylogeny/data/refseqs/rCRS.rCRS.fasta

	# [2/2] Prevent Bash from reading past this point once script is done
	exit $?
}
