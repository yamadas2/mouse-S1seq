#!/bin/bash

# Requirement

# Software
#conda create -y --no-default-packages \
#  -n mouse-S1seq \
#  -c bioconda \
#  FASTQC \
#  trim-galore \
#  bowtie2 \
#  samtools \
#  picard \
#  bedtools \
#  ucsc-bedgraphtobigwig \
#  ucsc-bigwigtobedgraph \
#  ucsc-wigtobigwig \
#  ucsc-bigwigtowig \
#  deeptools \
#  igvtools \
#  weblogo \
#  sra-tools \
#  hisat2 \
#  snpsplit
#  macs2
# Activate a conda environment
conda activate mouse-S1seq


# Files
# [GENOME].chrom.sizes files
# Bowtie2 indexes
# Peak bed files to create tag count matrixes around peaks with deepTools
# all_[GENOME].txt.gz SNP files for allele-specific allignment with SNPsplit


# Load environment variables
source ~/.bashrc


# Machine spec
CPU=${MAP_CPU:-4}
MEM=${MAP_MEM:-16}


# Temporary working directory
TMP=${MAP_TMP:-~/Desktop/TMP}
mkdir -p $TMP


# Directory for deep sequencing analyses
DIR_MAPPING=${MAP_DIR_MAPPING:-~/mapping}


# Data directory
DIR_FASTQ=$DIR_MAPPING/FASTQ/
DIR_MAP=$DIR_MAPPING/MAP/


# Chromosome sizes
DIR_CHROM_SIZES=${MAP_DIR_CHROM_SIZES:-$DIR_MAPPING/chrom.sizes}


# Bowtie2
export BOWTIE2_INDEXES=${MAP_BOWTIE2_INDEXES:-$TMP/indexes/}
if [ ! -d $BOWTIE2_INDEXES ]; then
  if [ ! -d $DIR_MAPPING/indexes ]; then
    echo "Bowtie2 indexes not found at $BOWTIE2_INDEXES or $DIR_MAPPING/indexes"
    exit
  else
    cp -r $DIR_MAPPING/indexes $BOWTIE2_INDEXES
  fi
fi


# Directories for deepTools and R
# peak bed files
DIR_PEAKS=${MAP_DIR_PEAKS:-$DIR_MAPPING/peaks}
# bigWig files
DIR_BWS=${MAP_DIR_BWS:-$DIR_MAPPING/bw}
# matrix files
DIR_MATS=${MAP_DIR_MATS:-$DIR_MAPPING/matrix.Rdata}
# rpm files
DIR_RPMS=${MAP_DIR_RPMS:-$DIR_MAPPING/rpm}
# plot files
DIR_PLOTS=${MAP_DIR_PLOTS:-$DIR_MAPPING/plots}


# Directory of SNP files for Allele-specific allignment with SNPsplit
DIR_SNPSPLIT=${MAP_DIR_SNPSPLIT:-$DIR_MAPPING/SNPsplit}


# Directory example
# mkdir -p $DIR_MAPPING \
#          $DIR_FASTQ \
#          $DIR_MAP \
#          $DIR_CHROM_SIZES \
#          $DIR_PEAKS/mm10/primary \
#          $DIR_PEAKS/hg19/primary \
#          $DIR_BWS \
#          $DIR_MATS \
#          $DIR_RPMS \
#          $DIR_PLOTS \
#          $DIR_SNPSPLIT
#          $DIR_MAPPING/indexes \
#          $DIR_MAPPING/scripts \


