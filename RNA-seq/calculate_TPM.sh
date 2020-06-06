
# テストデータ「count.txt」のダウンロード
wget https://gist.githubusercontent.com/akikuno/73882ef780dfafa5473a99cef76800fc/raw/c06b0f55652a65041ec58bce455397aa7388d314/count.txt

calc_rpk(){
    set /dev/stdin
    cat "${1}" |
    sed 1d |
    awk '{rpk=""
        for(i=3;i<=NF;i++){rpk=rpk" "$i/$2*1000}
        print $1, rpk}'
}

calc_scale(){
    set /dev/stdin
    cat "${1}" |
    tee temp_$$ |
    awk '{for(i=2;i<=NF;i++) sum[i]+=$i}
        END{for(key in sum) printf " "sum[key]}' |
    sed "s/^/sum /g" |
    awk 'NR==1 && FNR==1{ for(i=2;i<=NF;i++){denom[i]=$i}; next} {
        tpm=""
        for(i=2;i<=NF;i++){
            tpm=tpm","($i/denom[i] * 1000000)
            }
        print $1, tpm
        }' - temp_$$ |
    sed "s/ //g" |
    cat
    rm temp_$$
}

cat count.txt |
    calc_rpk |
    calc_scale |
cat > count_TPM.csv

# TPM check
cat count_TPM.csv |
    awk -F "," '{for(i=2;i<=NF;i++){sum[i]+=$i}}
    END{for(key in sum) print sum[key]}'

# RPKM
cat count.txt |
    calc_scale |
    calc_rpk |
cat > count_RPKM.csv
