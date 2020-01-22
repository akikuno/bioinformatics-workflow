# NGS-workflow

The latest [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installation is highly recommended.

```
conda update -y -n base conda
conda create -y -n ngs_analysis python=3.7 anaconda wget \
    r-base r-essentials r-tidyverse
conda install -y -n ngs_analysis -c bioconda \
    fastqc fastp samtools bedtools \
    star rsem kallisto r-sleuth bioconductor-deseq2 \
    bowtie2 macs2 deeptools homer \
    bioconductor-chippeakanno bioconductor-diffbind
conda activate ngs_analysis
git clone https://github.com/akikuno/ngs-workflow.git
```
