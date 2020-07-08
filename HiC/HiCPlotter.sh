#!/bin/sh

###############################################################################
#! Prepare software
###############################################################################
git clone https://github.com/kcakdemir/HiCPlotter.git

cat << 'EOF' > transpose.awk
{for (i=1; i<=NF; i++)  {a[NR,i] = $i};
} NF>p { p = NF }
END {for(j=1; j<=p; j++) {
        str=a[1,j];
        for(i=2; i<=NR; i++){str=str" "a[i,j]}
        print str
    }
}
EOF

###############################################################################
#! Download HiC data (size: 25.3GB)
###############################################################################

mkdir -p GSE72697_RAW
wget -P GSE72697_RAW -c "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72697/suppl/GSE72697_RAW.tar"
tar GSE72697_RAW -xf GSE72697_RAW.tar
mv GSM* GSE72697_RAW
find ./ -name "*matrix_cis*tar.gz" | xargs -I @ tar -xf @


###############################################################################
#! Formatting HiC matrix
###############################################################################

mkdir -p temp/

for input in $(find . -type d -name *iced-snpMasked* | grep txt | grep 40kb); do
    output=$(echo ${input} | sed "s/__/\t/g" | cut -f 4 | cut -d "/" -f 1)
    echo "${output} is now processing..."
    # Cast
    find . -type f -name "*chrX-cast.matrix.gz" |
        grep ${input} |
        xargs -I @ gzip -dc @ |
        awk -f transpose.awk |
    cat > temp/"${output}"_cast.matrix &
    # 129S1
    find . -type f -name "*chrX-129S1.matrix.gz" |
        grep ${input} |
        xargs -I @ gzip -dc @ |
        awk -f transpose.awk |
    cat > temp/"${output}"_129s1.matrix &
    wait 1>/dev/null 2>/dev/null
done

###############################################################################
#! Detect genome position in HiC
###############################################################################

# : > temp/tmp.bed
# printf "chrX\t4720001\t8320001\tInterest\n" >> temp/tmp.bed
# printf "chrX\t35758241\t35809632\tLAMP2\n" >> temp/tmp.bed
# printf "chrX\t103030677\t103170536\tAtrx\n" >> temp/tmp.bed
# cat temp/tmp.bed | sort -k 1,1 -k 2,2n > temp/genome_position.bed
# rm temp/tmp.bed

: > temp/tmp.bed
printf "chrX\t6500001\t8000001\tTarget\n" >> temp/tmp.bed
printf "chrX\t35000001\t36500001\tLAMP2\n" >> temp/tmp.bed
printf "chrX\t102500000\t104000000\tAtrx\n" >> temp/tmp.bed
cat temp/tmp.bed | sort -k 1,1 -k 2,2n > temp/genome_position.bed
rm temp/tmp.bed

head -n 1 temp/ESC_129s1.matrix |
    sed "s/ /\n/g" |
    awk '{print $0, NR}' |
    sed 1d |
    cut -d ":" -f 2- |
    sed "s/^/chrX /g" |
    sed "s/-/ /g" |
    sed "s/ /\t/g" |
    sort -k 1,1 -k 2,2n |
cat > temp/matrix.bed

bedtools intersect -wb -b temp/genome_position.bed -a temp/matrix.bed |
awk '{print $8,$4,$4}' |
awk '{if(min[$1] == "") min[$1]="inf"
if(min[$1]>$2) min[$1]=$2
if(max[$1]<$3) max[$1]=$3}
END{for(key in min){
    print key, min[key], max[key]
}}' |
cat > temp/target_loci.bed

###############################################################################
#! plot by HiCPlotter
###############################################################################

mkdir -p plot

#==============================================================================
#? All ChrX
#==============================================================================

for input in $(find . -type d -name *iced-snpMasked* | grep txt | grep 40kb | grep -v deltaFT)
do
    output=$(echo ${input} | sed "s/__/\t/g" | cut -f 4 | cut -d "/" -f 1)
    echo "${output} $1 is now processing..."

    python2 HiCPlotter/HiCPlotter.py \
        -f temp/"${output}"_cast.matrix temp/"${output}"_129s1.matrix \
        -n "Cast(Xa)" "129s1(Xi)" -chr chrX -o plot/${output}_ALL \
        -r 40000 -c 1 &
done 2>/dev/null
wait 2>/dev/null

#==============================================================================
#? Plot selected gene loci
#==============================================================================

# cat << EOF |
# Target 164 201
# Lamp2 886 906
# Atrx 2569 2589
# EOF
cat temp/target_loci.bed |
while read -r line; do
    set ${line}
    for input in $(find . -type d -name *iced-snpMasked* | grep txt | grep 40kb | grep -v deltaFT)
    do
        output=$(echo ${input} | sed "s/__/\t/g" | cut -f 4 | cut -d "/" -f 1)
        echo "${output} $1 is now processing..."

        python2 HiCPlotter/HiCPlotter.py \
            -f temp/"${output}"_cast.matrix temp/"${output}"_129s1.matrix \
            -n "Cast(Xa)" "129s1(Xi)" -chr chrX -o plot/${output}_"${1}"_"${2}"_"${3}"_LOG2 \
            -r 40000 -c 1 -s $2 -e $3 &

        python2 HiCPlotter/HiCPlotter.py \
            -f temp/"${output}"_cast.matrix temp/"${output}"_129s1.matrix \
            -n "Cast(Xa)" "129s1(Xi)" -chr chrX -o plot/${output}_"${1}"_"${2}"_"${3}"_TAD \
            -r 40000 -c 1 -s $2 -e $3 -ptd 1 -pi 1 &
    done 2>/dev/null
    wait 2>/dev/null
done