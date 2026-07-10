#!/usr/bin/env bash

MY_PWD=$1

mkdir $MY_PWD"/data"
cd $MY_PWD"/data"

# Download the reference genome for S. cerevisiae release 113
wget -O genome_sacCer3.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

# Unzip genome file
gunzip genome_sacCer3.fa.gz

# Download annotation file for S. cerevisiae release 113
wget -O genes_sacCer3.gtf.gz https://ftp.ensembl.org/pub/release-113/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.113.gtf.gz
gunzip genes_sacCer3.gtf.gz

# Create bowtie2 index from the reference genome
mkdir Bowtie2Index
docker run -i -v $MY_PWD"/data":$MY_PWD"/data" uschwartz/core_docker:v1.0 \
bowtie2-build $MY_PWD"/data/genome_sacCer3.fa" $MY_PWD"/data/Bowtie2Index/genome"

# Get the nucDetective pipeline
cd ..
git clone https://github.com/uschwartz/nucDetective.git
