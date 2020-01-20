# Install required packages by `conda`

conda create -y -n ngs-analysis python=3.7 anaconda r-base r-essentials r-tidyverse 
conda install -y -n ngs-analysis -c bioconda fastqc fastp hisat2 bowtie2 macs2 multiqc \
bioconductor-chipqc bioconductor-deseq2
