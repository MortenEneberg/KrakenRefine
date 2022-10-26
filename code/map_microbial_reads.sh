#!/usr/bin/env bash

module purge

Bowtie2=$1
reads_taxID=$2
sam=$3
report=$4

module load $Bowtie2

echo $reads_taxID

bowtie2 -x data/index/kraken2_genomes \
	-U $reads_taxID \
	-f \
	-p 20 \
	--local \
	--sensitive-local \
	-k 2 \
	--mp 2,2 \
	--score-min G,20,8 \
	--met-file $report \
	--no-hd \
	--mm \
	-S $sam

module purge
