# Install

## Update conda

```bash
conda update -y -n base conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## Create environments

```bash
conda info -e | grep -q "ngs$" || conda create -y -n ngs
conda install -y -n ngs -c conda-forge r-base r-essentials r-biocmanager
conda install -y -n ngs -c bioconda homer
```

```bash
PATH_HOMER=$(
  type homer |
  awk '{print $NF}' |
  sed "s|/bin/|/share/"
)

$PATH_HOMER/configureHomer.pl -install mm10 hg38
```
