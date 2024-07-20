# Installation

## Install software

- Install [VS code](https://azure.microsoft.com/ja-jp/products/visual-studio-code/)

- Install [R](https://cran.rstudio.com/) and [Rstudio](https://rstudio.com/products/rstudio/download/#download)

- Install [IGV](https://igv.org/doc/desktop/#DownloadPage/)


## Create directories for mapping
Download the [mapping](/mapping) directoroy


## Instration using Terminal

### Change the default shell from zsh to bash
Type the following in the Terminal
```bash
chsh -s /bin/bash
```

### Set environment variables
If a script is run with 4 CPUs and 16G memory and the "mapping" directory is located at the home directory, add the following lines to ~/.bashrc
```bash
MAP_CPU=4
MAP_MEM=16
MAP_DIR_MAPPING=~/mapping
```
Alternatively, if ~/.bashrc does not exist, download environment.variables.txt from GitHub, open the file, modify the variables and rename and copy the file to the home directory
```bash
cp environment.variables.txt ~/.bashrc
```

### Install [conda](https://docs.anaconda.com/miniconda/)
- apple silicon Mac
```bash
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

- Intel Mac
```bash
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

- Linux
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

### Initialize conda (for the first time only)
```bash
~/miniconda3/bin/conda init bash
```

### Create a conda environment and install packages
- Intel Mac
```bash
conda create -y --no-default-packages \
  -n mouse-S1seq \
  -c bioconda \
  FASTQC \
  trim-galore \
  bowtie2 \
  samtools \
  picard \
  bedtools \
  ucsc-bedgraphtobigwig \
  ucsc-bigwigtobedgraph \
  ucsc-wigtobigwig \
  ucsc-bigwigtowig \
  deeptools \
  igvtools \
  weblogo \
  macs2
```

- WSL2 Ubuntu22
Install packages by conda and apt-get, and Unicode locale for picard
```bash
conda create -y --no-default-packages \
  -n mouse-S1seq \
  -c bioconda \
  python=3.6 \
  FASTQC \
  python=2.7 \
  picard=2.18.7 \
  ucsc-bedgraphtobigwig=366 \
  ucsc-bigwigtobedgraph=366 \
  ucsc-wigtobigwig=366 \
  ucsc-bigwigtowig=366 \
  deeptools \
  igvtools=2.16.2 \
  weblogo \
  macs2
sudo apt-get install trim-galore \
  bowtie2 \
  samtools \
  r-base \
  bedtools
sudo locale-gen en_US.UTF-8
```

### Save conda package version information
Write conda package version information for reproducibility
```bash
conda list -n mouse-S1seq > conda.env.mouse-S1seq.txt
conda env export -n mouse-S1seq > conda.env.mouse-S1seq.yml
```
