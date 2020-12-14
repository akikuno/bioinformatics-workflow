#!/bin/sh

# ##################################################
# Define threads
# ##################################################
# Linux and similar...
[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/2)}')
# FreeBSD and similar...
[ -z "$threads" ] && threads=$(getconf NPROCESSORS_ONLN | awk '{print int($0/2)}')
# Solaris and similar...
[ -z "$threads" ] && threads=$(ksh93 -c 'getconf NPROCESSORS_ONLN' | awk '{print int($0/2)}')
# Give up...
[ -z "$threads" ] && threads=1


# ##################################################
# STAR index
# ##################################################

mkdir -p human_genome/human_index/STAR

STAR \
    --runMode genomeGenerate \
    --genomeDir human_genome/human_index/STAR \
    --genomeFastaFiles human_genome/*.primary_assembly.fa \
    --sjdbGTFfile human_genome/*.gtf \
    --runThreadN "$threads"

# ##################################################
# STAR mapping
# ##################################################

R1=$(ls ./fastq/*R1*.gz)
R2=$(ls ./fastq/*R2*.gz)
num=$(find ./fastq/*R1*.gz -type f | awk '{print NR}')

mkdir -p bam

for i in $(echo $num) ; do
    fw=$(echo $R1 | cut -d " " -f $i)
    rv=$(echo $R2 | cut -d " " -f $i)
    out_f=$(echo "$fw" |
    sed -e "s#.*/#bam/#g" -e "s/_R1.*//g")
    #
    time STAR --runThreadN "$threads" \
        --genomeDir human_genome/human_index/STAR \
        --readFilesIn "$fw" "$rv" \
        --readFilesCommand gunzip -c \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix "$out_f"_
done

# ##################################################
# Sort and index
# ##################################################

for bam in ./bam/*bam ; do
    out_f=$(echo $bam | sed "s/\.bam/_sorted.bam/g")
    samtools sort -@ "$threads" $bam -o "$out_f"
    samtools index -@ "$threads" "$out_f"
    rm $bam
done
