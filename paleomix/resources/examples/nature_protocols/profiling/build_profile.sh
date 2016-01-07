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

OUTPUT_SAM=${OUTPUT_PREFIX}.bowtie2.out.sam
OUTPUT_SORTED_BAM=${OUTPUT_PREFIX}.sorted.bam
OUTPUT_RMDUP_BAM=${OUTPUT_PREFIX}.noduplicates.bam
OUTPUT_METAPHLAN_INPUT=${OUTPUT_PREFIX}.noduplicates
OUTPUT_METAPHLAN=${OUTPUT_PREFIX}.txt

if [ ! -e "${JAR_ROOT}/SortSam.jar" ];
then
    echo "Required JAR file is missing; please install:"
    echo "  - ${JAR_ROOT}/SortSam.jar"
    exit 1
elif [ ! -e "${METAPHLAN_ROOT}/metaphlan.py" ];
then
    echo "MetaPhlAn does not appear to be installed; please install at:"
    echo "  - ${METAPHLAN_ROOT}/metaphlan.py"
    exit 1
fi

for executable in bzcat bowtie2 java bam_rmdup_collapsed samtools;
do
    if ! which ${executable} > /dev/null;
    then
	echo "Required executable is missing ('${executable}'); please install."
	exit 1
    fi
done


echo
echo "Generating profile from $# files, saving to ${OUTPUT_METAPHLAN}"

if [ ! -e "${OUTPUT_METAPHLAN}" ];
then
    if [ ! -e "${OUTPUT_METAPHLAN_INPUT}" ];
    then
	if [ ! -e "${OUTPUT_RMDUP_BAM}" ];
	then
	    if [ ! -e "${OUTPUT_SORTED_BAM}" ];
	    then
		if [ ! -e "${OUTPUT_SAM}" ];
		then
		    bzcat $@ | bowtie2 -x "${BOWTIE_PREFIX}" -U - -S "${OUTPUT_SAM}.tmp" --no-unal --quiet
		    mv "${OUTPUT_SAM}.tmp" "${OUTPUT_SAM}"
		fi

		java -Xmx4g -jar ${JAR_ROOT}/SortSam.jar SO=coordinate VERBOSITY=WARNING QUIET=TRUE \
		    "I=${OUTPUT_SAM}" "O=${OUTPUT_SORTED_BAM}.tmp"
		mv "${OUTPUT_SORTED_BAM}.tmp" "${OUTPUT_SORTED_BAM}"
	    fi

	    bam_rmdup_collapsed --remove-duplicates < "${OUTPUT_SORTED_BAM}" > "${OUTPUT_RMDUP_BAM}.tmp"
	    mv "${OUTPUT_RMDUP_BAM}.tmp" "${OUTPUT_RMDUP_BAM}"
	fi

	samtools view ${OUTPUT_RMDUP_BAM} | awk '{print $1 "\t" $3}' > "${OUTPUT_METAPHLAN_INPUT}.tmp";
	mv "${OUTPUT_METAPHLAN_INPUT}.tmp" "${OUTPUT_METAPHLAN_INPUT}"
    fi

    ${METAPHLAN_ROOT}/metaphlan.py "${OUTPUT_METAPHLAN_INPUT}" > "${OUTPUT_METAPHLAN}.tmp"
    mv "${OUTPUT_METAPHLAN}.tmp" "${OUTPUT_METAPHLAN}"
fi

echo "Done: Profile written to ${OUTPUT_METAPHLAN}"
