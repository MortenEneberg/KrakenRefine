#!/usr/bin/env bash

module purge

Kraken2=$1
kraken2db=$2
confidence=$3
Outdir=$4
reads=$5
output_report=$6
output_kraken=$7

mkdir $4

