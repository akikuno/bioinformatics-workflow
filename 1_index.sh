#!/bin/sh

# Define threads
# Linux and similar...
[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/2)}')
# FreeBSD and similar...
[ -z "$threads" ] && threads=$(getconf NPROCESSORS_ONLN | awk '{print int($0/2)}')
# Solaris and similar...
[ -z "$threads" ] && threads=$(ksh93 -c 'getconf NPROCESSORS_ONLN' | awk '{print int($0/2)}')
# Give up...
[ -z "$threads" ] && threads=1

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# STAR, RSEM, kallisto, bowtie2 index
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mkdir -p mouse_index/STAR mouse_index/RSEM mouse_index/kallisto mouse_index/bowtie2

STAR \
--runMode genomeGenerate \
--genomeDir mouse_index \
--genomeFastaFiles mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile mouse_genome/Mus_musculus.GRCm38.99.gtf \
--runThreadN "$threads"

rsem-prepare-reference \
--gtf mouse_genome/Mus_musculus.GRCm38.99.gtf \
--num-threads "$threads" \
mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
mouse_index/RSEM/index

kallisto index -i mouse_index/kallisto \
mouse_genome/Mus_musculus.GRCm38.cdna.all.fa

bowtie2-build \
mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
mouse_index/bowtie2/index \
--threads "$threads"