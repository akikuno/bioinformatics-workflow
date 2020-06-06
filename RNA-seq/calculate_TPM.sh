
# テストデータ「count.txt」のダウンロード
wget -O - https://gist.githubusercontent.com/akikuno/73882ef780dfafa5473a99cef76800fc/raw/c06b0f55652a65041ec58bce455397aa7388d314/count.txt |
cat > count.txt

calc_rpk(){
    set /dev/stdin
    cat "${1}" |
    awk 'NR==1{$2=""; print $0}
        NR>1{rpk=""
        for(i=3;i<=NF;i++){rpk=rpk" "$i/$2*1000}
        print $1, rpk}'
    sed "s/  */\t/g" |
}

calc_scale(){
    set /dev/stdin
    cat "${1}" |
    tee temp_$$ |
    awk 'NR==1
        NR>1{for(i=2;i<=NF;i++) sum[i]+=$i}
        END{for(key in sum) printf " "sum[key]}' |
    # awk '{print $0, FNR, NR , FILENAME}' - temp_$$
    awk 'NR==1 && FNR==1
        NR==2 && FNR==2 { for(i=1;i<=NF;i++) denom[i]=$i}
        NR>3{
        tpm=""
        for(i=2;i<=NF;i++){
            tpm=tpm" "($(i)/denom[i-1] * 1000000)
            }
        print $1, tpm
        }
        ' - temp_$$ |
    sed "s/  */\t/g" |
    cat
    rm temp_$$
}

cat count.txt |
    calc_rpk |
    calc_scale |
cat > count_TPM.txt

# TPM check
cat count_TPM.txt |
    awk '{for(i=2;i<=NF;i++){sum[i]+=$i}}
    END{for(key in sum) print sum[key]}'

# RPKM
cat count.txt |
    calc_scale |
    calc_rpk |
cat > count_RPKM.txt
