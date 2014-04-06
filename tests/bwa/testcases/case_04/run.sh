#!/bin/bash
set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

cd $(dirname $0)

BWA=$1
PREFIX=prefix.fasta

rm -f ${PREFIX}.*

$BWA index ${PREFIX} 2> run.log
$BWA aln ${PREFIX} reads.fasta > reads.fasta.fai 2> run.log

$BWA samse ${PREFIX} reads.fasta.fai reads.fasta 2> run.log > results.sam

set +o pipefail
samtools view -S results.sam 2>&1|grep error > run.log && exit 1 || exit 0
