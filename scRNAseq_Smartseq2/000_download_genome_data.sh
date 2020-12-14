#!/bin/sh

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Download human genome data: hg38
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

mkdir -p human_genome

# Whole-genome sequence for STAR mapping
wget -qO - \
ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz |
gzip -cd > human_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa &

# # cDNA for kallisto mapping
# wget -qP human_genome \
# ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz &&
# gzip -d mouse_genome/Mus_musculus.GRCm38.cdna.all.fa.gz &

# Gene information
wget -qO - \
ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz |
gzip -cd > human_genome/Homo_sapiens.GRCh38.100.gtf &

wait

rm wget-log*

exit 0
