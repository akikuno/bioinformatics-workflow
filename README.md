# NGS-workflow

The latest [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installation is highly recommended.

## Channel setup
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
## Create environment and install software
```
conda update -y -n base conda
conda create -y -n ngs_analysis python=3.7 anaconda wget \
    r-base r-essentials r-tidyverse r-ggpubr \
    fastqc fastp samtools bedtools \
    star rsem kallisto r-sleuth bioconductor-deseq2 \
    bowtie2 macs2 deeptools homer \
    bioconductor-chippeakanno bioconductor-diffbind
conda activate ngs_analysis
```
