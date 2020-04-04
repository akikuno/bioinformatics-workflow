#!/bin/sh

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Download mouse genome data: mm10, GRCm38
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

mkdir -p mouse_genome

# Whole-genome sequence for STAR mapping
wget -qP mouse_genome \
ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz &&
gzip -d mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz &

# cDNA for kallisto mapping
wget -qP mouse_genome \
ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz &&
gzip -d mouse_genome/Mus_musculus.GRCm38.cdna.all.fa.gz &

# Gene information
wget -qP mouse_genome \
ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz &&
gzip -d mouse_genome/Mus_musculus.GRCm38.99.gtf.gz &

wait

rm wget-log*

exit 0
