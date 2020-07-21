# NGS-workflow

The latest [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installation is highly recommended.

## 1. Channel setup
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## 2. Create environment and install software
```bash
if [ $(conda info -e | cut -d " " -f 1 | grep -c bfx$) -eq 0 ]
then
conda update -y -n base conda
conda create -y -n bfx python=3.7 anaconda wget \
    r-base r-essentials r-ggpubr \
    fastqc fastp qualimap multiqc \
    samtools bedtools deeptools \
    star rsem kallisto subread r-sleuth r-seurat \
    bioconductor-deseq2 \
    bwa bowtie bowtie2 macs2 homer \
    bioconductor-chippeakanno bioconductor-diffbind
fi
```

## 3. Activate the environment
```bash
conda activate bfx
```
You can skip step 1 and 2 from the second time on.

