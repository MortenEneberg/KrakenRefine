# README
KrakenRefine
=====
A Snakemake pipeline that takes **Kraken** (Kraken2) results and **Refine** them by aligning a subset of reads of particular interest with Bowtie2 to a set of reference genomes. The presence of a particular group of species (natively grouped at genus level, but this can be changed), is then evaluated based on coverage profiles, because many false positives origin from database genome contamination that lead to anomalies in coverage profiles.  

It was developed to reeva
pipeline to proces kraken2 data to remove false positives by mapping reads to reference genomes
This is a test