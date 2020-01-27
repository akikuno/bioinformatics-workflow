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

##----------------------------------------------
# Format FASTQ files
##----------------------------------------------

#= parameters ==============
dir="FASTQ DIR PATH"
#===========================

mkdir -p fastq
for subdir in $(ls "$dir") ; do
    cat "$dir"/"$subdir"/*R1* > fastq/"$subdir"_R1.fq.gz &
    cat "$dir"/"$subdir"/*R2* > fastq/"$subdir"_R2.fq.gz &
    wait 2>/dev/null 1>/dev/null
done

## FASTQC
mkdir -p ./fastqc_report
for file in ./fastq/*; do
    fastqc -t "$threads" "$file" \
    --nogroup -o ./fastqc_report
done

##----------------------------------------------
# Multiqc
##----------------------------------------------
multiqc .

##----------------------------------------------
# fastp trimming
##----------------------------------------------

R1=$(ls ./fastq/*R1.fq.gz)
R2=$(ls ./fastq/*R2.fq.gz)

num=$(find ./fastq/*R1.fq.gz -type f | awk '{print NR}')
#awk '{if (NR%2!=0) a++; print a }'

mkdir -p fastq_trim

for i in $(echo $num) ; do
    fw=$(echo $R1 | cut -d " " -f $i)
    rv=$(echo $R2 | cut -d " " -f $i)
    echo $fw
    
    out_fw=$(echo ${fw} | sed "s/fastq/fastq_trim/g")
    out_rv=$(echo ${rv} | sed "s/fastq/fastq_trim/g")
    report=$(echo ${out_fw} | sed "s/_R1.fq.gz//g")
    #
    time fastp -i ${fw} -I ${rv} \
    -o ${out_fw} -O ${out_rv} \
    --trim_front1 5 -trim_tail1 1 \
    -h ${report}.html -j ${report}.json \
    --thread ${threads}
done

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
