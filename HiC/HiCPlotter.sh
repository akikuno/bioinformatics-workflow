#!/bin/sh

# ===========================================================
# Prepare software
# ===========================================================
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

# ===========================================================
# Download HiC data (size: 25.3GB)
# ===========================================================
mkdir -p GSE72697_RAW
wget -P GSE72697_RAW -c "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72697/suppl/GSE72697_RAW.tar"
tar GSE72697_RAW -xf GSE72697_RAW.tar
mv GSM* GSE72697_RAW
find ./ -name "*matrix_cis*tar.gz" | xargs -I @ tar -xf @

# ===========================================================
# Visualize HiC data on X chromosome
# ===========================================================

for input in $(find . -type d -name *iced-snpMasked* | grep txt | grep 40kb); do
    output=$(echo ${input} | sed "s/__/\t/g" | cut -f 4 | cut -d "/" -f 1)
    echo "${output} is now processing..."
    # -------------------------------------------
    # Formatting HiC matrix
    # -------------------------------------------
    # Cast
    find . -type f -name "*chrX-cast.matrix.gz" |
    grep ${input} |
    xargs -I @ gzip -dc @ |
    awk -f transpose.awk \
    > tmp_cast.matrix &
    # 129S1
    find . -type f -name "*chrX-129S1.matrix.gz" |
    grep ${input} |
    xargs -I @ gzip -dc @ |
    awk -f transpose.awk \
    > tmp_129s1.matrix &
    wait 1>/dev/null 2>/dev/null
    # -------------------------------------------
    # plot by HiCPlotter
    # -------------------------------------------
    # mm9: chrX:4,720,001-8,320,001
    python2 HiCPlotter/HiCPlotter.py \
    -f tmp_cast.matrix tmp_129s1.matrix \
    -n "Cast(Xa)" "129s1(Xi)" -chr chrX -o ${output}_LOG2 \
    -r 40000 -c 1 -s 120 -e 210 &
    # mm9: chrX:4,720,001-8,320,001
    python2 HiCPlotter/HiCPlotter.py \
    -f tmp_cast.matrix tmp_129s1.matrix \
    -n "Cast(Xa)" "129s1(Xi)" -chr chrX -o ${output}_TAD \
    -r 40000 -c 1 -s 120 -e 210 -ptd 1 -pi 1 &
    # mm9: chrX:4,720,001-8,320,001
    python2 HiCPlotter/HiCPlotter.py \
    -f tmp_cast.matrix tmp_129s1.matrix \
    -n "Cast(Xa)" "129s1(Xi)" -chr chrX -o ${output}_ALL \
    -r 40000 -c 1 &
    wait 1>/dev/null 2>/dev/null
    rm tmp_*
done

