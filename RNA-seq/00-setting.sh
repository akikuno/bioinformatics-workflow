#!/bin/bash

# Make directories
mkdir -p data/bam data/bigwig data/fastq/ data/mouse_genome data/mouse_index reports/fastqc

echo "*" >data/.gitignore

# Install required software

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -y -n ngs
conda install -y -n ngs wget fastqc sra-tools subread deeptools r-base r-essentials
