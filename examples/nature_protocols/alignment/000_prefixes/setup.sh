#!/bin/bash

set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

cd $(dirname $0)

if [ ! -e "Pi_nucl.fasta" ];
then
    echo "Fetching P. infestans T30-4 nuclear genome"
    echo "  See http://protists.ensembl.org/Phytophthora_infestans/Info/Index"
    curl "ftp://ftp.ensemblgenomes.org/pub/protists/release-20/fasta/phytophthora_infestans/dna/Phytophthora_infestans.ASM14294v1.20.dna.toplevel.fa.gz" \
	-C - -o "Pi_nucl.fasta.gz"
    gunzip "Pi_nucl.fasta.gz"
else
    echo "Pi_mito.fasta already fetched; skipping ..."
fi

echo
if [ ! -e "Pi_mito.fasta" ];
then
    echo "Fetching P. infestans mitochondrial genome"
    echo "  See http://www.ncbi.nlm.nih.gov/nuccore/AY894835.1"
    curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AY894835.1&rettype=fasta&retmode=text" \
	-C - -o "Pi_mito.fasta"
else
    echo "Pi_mito.fasta already fetched; skipping ..."
fi


echo
if [ ! -e "St_nucl.fasta" ];
then
    echo "Fetching S. tuberosum nuclear genome"
    echo "  See http://solanaceae.plantbiology.msu.edu/"
    curl "http://solanaceae.plantbiology.msu.edu/data/PGSC_DM_v3_superscaffolds.fasta.zip" \
	-C - -o "St_nucl.fasta.zip"
    unzip "St_nucl.fasta.zip"
    rm "St_nucl.fasta.zip"
    mv -v PGSC_DM_v3_superscaffolds.fasta St_nucl.fasta
else
    echo "St_nucl.fasta already fetched; skipping ..."
fi

echo
if [ ! -e "St_mito.fasta" ];
then
    echo "Fetching S. tuberosum mitochondrial genome"
    echo "  See http://solanaceae.plantbiology.msu.edu/"
    curl "http://solanaceae.plantbiology.msu.edu/data/S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta.zip" \
	-o "St_mito.fasta.zip"
    unzip "St_mito.fasta.zip"
    rm "St_mito.fasta.zip"
    mv -v S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta St_mito.fasta
else
    echo "St_mito.fasta already fetched; skipping ..."
fi

echo
if [ ! -e "St_chlor.fasta" ];
then
    echo "Fetching S. tuberosum chloroplast genome"
    echo "  See http://solanaceae.plantbiology.msu.edu/"
    curl "http://solanaceae.plantbiology.msu.edu/data/S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta.zip" \
	-o "St_chlor.fasta.zip"
    unzip "St_chlor.fasta.zip"
    rm "St_chlor.fasta.zip"
    mv -v S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta St_chlor.fasta
else
    echo "St_chlor.fasta already fetched; skipping ..."
fi
