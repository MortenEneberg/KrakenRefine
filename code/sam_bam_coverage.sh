#!/usr/bin/env bash

module purge

Samtools=$1
Bedtools=$2
sam_header=$3
bam=$4
bed=$5

module load $Samtools

samtools view -S -b $sam_header | sort > $bam

module purge

module load $Bedtools

bedtools genomecov -ibam $bam > $bed

module purge