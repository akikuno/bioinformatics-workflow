# NGS-workflow

The latest [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installation is required.

## 1. Channel setup

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## 2. Create environment and install software

```bash
if [ $(conda info -e | cut -d " " -f 1 | grep -c "ngs$") -eq 0 ] ; then
conda update -y -n base conda
conda create -y -n ngs python=3.9 wget \
  sra-tools fastqc fastp qualimap multiqc \
  samtools bedtools deeptools \
  star subread bwa macs2 homer \
  r-base r-essentials r-ggpubr \
  r-seurat bioconductor-deseq2 \
  bioconductor-chippeakanno bioconductor-diffbind
fi
```

## 3. Activate the environment
```bash
conda activate ngs
```
You can skip step 1 and 2 from the second time on.

