#!/bin/bash

set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

cd $(dirname $0)

for ENA in `ls */000_ENA | sort`;
do
    SAMPLE=$(dirname ${ENA})
    echo "Fetching FASTQ reads for sample '${SAMPLE}' ..."

    tail -n +2 ${ENA} | cut -f10 | tr ";" "\n" |
    while read URL;
    do
	FNAME=$(echo $URL | sed -e's#.*/##')
	echo "  - $FNAME"

	if [ ! -e "${SAMPLE}/${FNAME}.DONE" ];
	then
	    PREFIX=${FNAME/.fastq.gz/}
	    # Split into chunks of 10M reads
	    curl "${URL}" | gunzip | split -l 40000000 -a 3 - "${SAMPLE}/${PREFIX}_"

	    ls ${SAMPLE}/${PREFIX}_* |
	    while read FNAME;
	    do
		mv ${FNAME} ${FNAME}.fastq
	    done

	    ls ${SAMPLE}/${PREFIX}_*.fastq | xargs -n 1 -P 8 gzip --verbose

	    touch "${SAMPLE}/${FNAME}.DONE";
	fi
    done
done


