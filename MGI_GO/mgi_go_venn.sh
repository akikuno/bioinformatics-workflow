#!/bin/sh

# http://www.informatics.jax.org/downloads/reports/index.html#go

# GO Evidence Codeについて:
# http://geneontology.org/docs/guide-go-evidence-codes/

mkdir -p mgi_go/data mgi_go/gene_list

# =================================================================
# MGI RNA binding/splicing 遺伝子リスト
# =================================================================

wget -P mgi_go/data http://www.informatics.jax.org/downloads/reports/gene_association.mgi.gz

rna_binding="GO:0003723"
rna_splicing="GO:0008380"

# -------------------------------------------
# Gene list associated with RNA binding
# -------------------------------------------
gzip -dc mgi_go/data/gene_association.mgi.gz |
    awk -v go="$rna_binding" '$4==go' | 
    cut -f 3 |
    sort -u |
cat - > mgi_go/gene_list/rna_binding.txt

# -------------------------------------------
# Gene list associated with RNA splicing
# -------------------------------------------
gzip -dc mgi_go/data/gene_association.mgi.gz |
    awk -v go="$rna_splicing" '$4==go' |
    cut -f 3 |
    sort -u |
cat - > mgi_go/gene_list/rna_splicing.txt

wc -l mgi_go/gene_list/rna_*

# =================================================================
# MGI embryonic lethality 遺伝子リスト
# =================================================================

wget -P mgi_go/data http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt

# -------------------------------------------
# Gene list associated with emb lethality E9.5-E14.5
# -------------------------------------------
MP="MP:0006207 MP:0011098 MP:0011108"

echo mgi_go/data/MGI_PhenoGenoMP.rpt |
    awk -v mp="$MP" \
    '{gsub("MP","-e MP",mp)
    print "grep",mp, $0}' |
    sh - |
    cut -d "<" -f 1 |
    cut -d "/" -f 1 |
    cut -f 1 |
    sort -u |
cat - > mgi_go/gene_list/emb_lethal_e14.txt

# -------------------------------------------
# Gene list associated with emb lethality E14.5-E18.5
# -------------------------------------------

MP="MP:0006208 MP:0011099 MP:0011109"

echo mgi_go/data/MGI_PhenoGenoMP.rpt |
    awk -v mp="$MP" \
    '{gsub("MP","-e MP",mp)
    print "grep",mp, $0}' |
    sh - |
    cut -d "<" -f 1 |
    cut -d "/" -f 1 |
    cut -f 1 |
    sort -u |
cat - > mgi_go/gene_list/emb_lethal_e18.txt

wc -l mgi_go/gene_list/emb_lethal*

# ちなみにe14とe18での重複遺伝子数は**228遺伝子**
join mgi_go/gene_list/emb_lethal_e14.txt mgi_go/gene_list/emb_lethal_e18.txt |
wc -l

# =================================================================
# venn diagram
# =================================================================
# RNA bindingとe14
join mgi_go/gene_list/rna_binding.txt mgi_go/gene_list/emb_lethal_e14.txt |
cat - > mgi_go/gene_list/venn_rna_binding_e14.txt

# RNA splicingとe14
join mgi_go/gene_list/rna_splicing.txt mgi_go/gene_list/emb_lethal_e14.txt |
cat - > mgi_go/gene_list/venn_rna_splicing_e14.txt

# RNA bindingとe18
join mgi_go/gene_list/rna_binding.txt mgi_go/gene_list/emb_lethal_e18.txt |
cat - > mgi_go/gene_list/venn_rna_binding_e18.txt

# RNA splicingとe18
join mgi_go/gene_list/rna_splicing.txt mgi_go/gene_list/emb_lethal_e18.txt |
cat - > mgi_go/gene_list/venn_rna_splicing_e18.txt

wc -l mgi_go/gene_list/venn_*

# =================================================================
# レポート用にひとつのファイルに集計
# =================================================================

true > mgi_go/resut.csv
for input in mgi_go/gene_list/venn*; do
    go=$(echo "$input" | sed -e "s/.*venn_//g" -e "s/_e.*.txt//g")
    day=$(echo "$input" | sed -e "s/.*_e/e/g" -e "s/.txt//g")
    #
    cat $input |
        sed "s/^/$go,$day,/g" |
    cat - >> mgi_go/resut.csv
done