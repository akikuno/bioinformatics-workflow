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

# # ##################################################
# # Quantification by RSEM
# # ##################################################

# # Index
# mkdir -p mouse_index/RSEM/index
# rsem-prepare-reference \
# --star \
# --gtf mouse_genome/Mus_musculus.GRCm38.99.gtf \
# --num-threads "$threads" \
# mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
# mouse_index/RSEM/index

# # Quantification
# fastqs="fastq fastq_trim"
# bams="bam bam_trim"
# mkdir -p ${bams}

# for ijk in 1 2; do
#     bam=$(echo ${bams} | cut -d " " -f ${ijk})
#     fastq=$(echo ${fastqs} | cut -d " " -f ${ijk})
#     #
#     R1=$(ls ./${fastq}/*R1*.gz)
#     R2=$(ls ./${fastq}/*R2*.gz)
#     num=$(find ./${fastq}/*R1*.gz -type f | awk '{print NR}')
#     #
#     for i in $(echo ${num}) ; do
#         fw=$(echo ${R1} | cut -d " " -f ${i})
#         rv=$(echo ${R2} | cut -d " " -f ${i})
#         out_f=$(echo "${fw}" |
#         sed -e "s/${fastq}/${bam}/g" -e "s/_R1*.gz//g")
#         echo "${out_f} is now processing..."
#         #
#         { gzip -dc ${fw} > /tmp/R1* & } 2>/dev/null
#         { gzip -dc ${rv} > /tmp/R2* & } 2>/dev/null
#         wait 2>/dev/null
#         time rsem-calculate-expression \
#         -p ${threads} \
#         --paired-end \
#         --star \
#         --output-genome-bam \
#         /tmp/R1* /tmp/R2* \
#         mouse_index/RSEM/index \
#         ${out_f}
#     done
# done
# multiqc .
