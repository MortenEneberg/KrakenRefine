#!/usr/bin/env bash

module purge

Kraken2 = $1
kraken2db = $2
confidence = $3
Outdir = $4

mkdir $Outdir

module load $Kraken2
