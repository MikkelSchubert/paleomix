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

if [ $# -lt 2 ];
then
    echo "Usage: $0 <output-prefix> <input_1.bz2> [<input_2.bz2>, ...]"
    exit 1
fi

OUTPUT_PREFIX=$1
# Remove the prefix parameter from the list of parameters
shift 1

JAR_ROOT=~/install/jar_root
METAPHLAN_ROOT=~/install/metaphlan
BOWTIE_PREFIX=${METAPHLAN_ROOT}/bowtie2db/mpa

echo "Generateing profile from $# files, saving to ${OUTPUT_PREFIX}:"
echo "  - Mapping $# files against MetaPhlan database at ${BOWTIE_PREFIX} ..."
bzcat $@ | bowtie2 -x "${BOWTIE_PREFIX}" -U - -S ${OUTPUT_PREFIX}.bowtie2.out.sam --no-unal

echo "  - Sorting reads ..."
java -Xmx4g -jar ${JAR_ROOT}/SortSam.jar I=${OUTPUT_PREFIX}.bowtie2.out.sam O=${OUTPUT_PREFIX}.sorted.bam SO=coordinate

echo "  - Removing PCR duplicates ..."
bam_rmdup_collapsed --remove-duplicates < ${OUTPUT_PREFIX}.sorted.bam > ${OUTPUT_PREFIX}.noduplicates.bam

echo "  - Collecting names of reads / markers ..."
samtools view ${OUTPUT_PREFIX}.noduplicates.bam | awk '{print $1 "\t" $3}' > ${OUTPUT_PREFIX}.noduplicates

echo "  - Building MetaPhlAn profile ..."
${METAPHLAN_ROOT}/metaphlan.py ${OUTPUT_PREFIX}.noduplicates > ${OUTPUT_PREFIX}.txt

echo "Done: Wrote profile to ${OUTPUT_PREFIX}.txt"
