#!/bin/bash
set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

cd $(dirname $0)

BWA=$1
PREFIX=prefix.fasta

rm -f ${PREFIX}.*

$BWA index ${PREFIX} 2> ${PREFIX}.log
$BWA aln ${PREFIX} reads1.fasta > reads1.fasta.fai 2> reads1.fasta.log
$BWA aln ${PREFIX} reads2.fasta > reads2.fasta.fai 2> reads2.fasta.log

set +o pipefail # Fail is a command in a chain of pipes fails
$BWA sampe ${PREFIX} reads1.fasta.fai reads2.fasta.fai reads1.fasta reads2.fasta 2>&1 | grep "Assertion" > run.log && exit 13 || exit 0

#../../check_sam.py --flags-unset 4 --position -1 results.sam