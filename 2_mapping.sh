#!/bin/sh

##----------------------------------------------
# STAR mapping
##----------------------------------------------

#= parameters ==============
star_index="PATH"
#===========================

R1=$(ls ./fastq/*R1.fq.gz)
R2=$(ls ./fastq/*R2.fq.gz)

num=$(find ./fastq/*R1.fq.gz -type f | awk '{print NR}')
#awk '{if (NR%2!=0) a++; print a }'

mkdir -p bam

for i in $(echo $num) ; do
    fw=$(echo $R1 | cut -d " " -f $i)
    rv=$(echo $R2 | cut -d " " -f $i)
    out_f=$(echo "$fw" |
    sed -e "s/fastq/bam/g" -e "s/_R1.fq.gz//g")
    
    time STAR --runThreadN "$threads" \
    --genomeDir "$star_index" \
    --readFilesIn "$fw" "$rv" \
    --readFilesCommand gunzip -c \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "$out_f"_
done

##----------------------------------------------
# Sort and index
##----------------------------------------------

for bam in ./bam/*bam ; do
    out_f=$(echo $bam | sed "s/\.bam/_sorted.bam/g")
    samtools sort -@ "$threads" $bam -o "$out_f"
    samtools index -@ "$threads" "$out_f"
    rm $bam
done

##----------------------------------------------
# BigWig files to visualize by IGV
##----------------------------------------------

mkdir -p bw

for bam in ./bam/*bam ; do
    out_f=$(echo "$bam" | sed -e "s/bam/bw/g" -e "s/sorted.bw/sorted_bin5_cpm.bw/g")
    time bamCoverage -b $bam -o "$out_f" \
    -p "$threads" \
    --binSize 5 \
    --normalizeUsing CPM
done
