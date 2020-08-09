#!/bin/sh

###############################################################################
#! Prerequisits
###############################################################################
# Python 2.7
# bedtools

###############################################################################
#! Download HiCPlotter
###############################################################################
[ -d HiCPlotter ] || git clone https://github.com/kcakdemir/HiCPlotter.git

###############################################################################
#! Download HiC data (size: 25.3GB)
###############################################################################

if [ -f GSE72697_RAW.tar ]; then
    mkdir -p GSE72697_RAW
    wget -P GSE72697_RAW -c "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72697/suppl/GSE72697_RAW.tar"
    tar GSE72697_RAW -xf GSE72697_RAW.tar
    mv GSM* GSE72697_RAW
    find ./ -name "*matrix_cis*tar.gz" | xargs -I @ tar -xf @
fi

###############################################################################
#! Formatting HiC matrix
###############################################################################

mkdir -p temp/

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

cat << EOF |
chrX 6000001 8000001 Target_02
chrX 4000000 24000001 Target_20
EOF
    sed "s/ /\t/g" |
    sort -k 1,1 -k 2,2n |
cat > temp/genome_position.bed

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

for output in $(find temp/*matrix | cut -d "_" -f 1 | cut -d "-" -f 1 | sed "s#.*/##g" | sort -u)
do
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

cat temp/target_loci.bed |
while read -r line; do
    set ${line}
    for output in $(find temp/*matrix | cut -d "_" -f 1 | cut -d "-" -f 1 | sed "s#.*/##g" | sort -u)
    do
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