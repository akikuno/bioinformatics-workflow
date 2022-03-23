#!/bin/sh

###############################################################################
# Download mouse genome data: GRCm38 (Ensembl release 102)
###############################################################################

mkdir -p data/mouse_genome

release=102

# Whole-genome sequence for mapping
wget -qP data/mouse_genome \
    ftp://ftp.ensembl.org/pub/release-"$release"/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz &&
    gzip -d data/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz &

# cDNA for kallisto mapping
wget -qP data/mouse_genome \
    ftp://ftp.ensembl.org/pub/release-"$release"/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz &&
    gzip -d data/mouse_genome/Mus_musculus.GRCm38.cdna.all.fa.gz &

# Gene information
wget -qP data/mouse_genome \
    ftp://ftp.ensembl.org/pub/release-"$release"/gtf/mus_musculus/Mus_musculus.GRCm38."$release".gtf.gz &&
    gzip -d data/mouse_genome/Mus_musculus.GRCm38."$release".gtf.gz &

wait

rm wget-log*
