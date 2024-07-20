#!/bin/bash

SDIR="$( dirname "$( cd "$( dirname "$0" )" && pwd )" )"


if [ "$#" -lt 3 ]; then
  echo "Usage mapping.s1seq.sh PROJECT NAME GENOME ADAPTOR_SWAP"
  exit
elif [ "$#" -gt 3 -a $4 != "TRUE" -a $4 != "FALSE" ]; then
  echo "Assign TRUE or FALSE to ADAPTOR_SWAP"
  exit
fi


PROJECT=$1
NAME=$2
GENOME=$3
ADAPTOR_SWAP=${4:-FALSE}


# Environment setting
source "$SDIR"/lib/mapping_variables.20240711.sh


set -uex


echo "-------------------------------------"
echo "Project      : $PROJECT"
echo "Sample       : $NAME"
echo "Genome       : $GENOME"
echo "Adaptor swap : $ADAPTOR_SWAP"
echo "CPU core     : $CPU"
echo "Memory size  : $MEM"
echo "Tmp folder   : $TMP"
echo
: $PROJECT $NAME $GENOME $ADAPTOR_SWAP $CPU $MEM $TMP


# Quality check with FastQC
DIR_FASTQC=$DIR_MAP/$PROJECT/fastqc
mkdir -p $DIR_FASTQC
fastqc $DIR_FASTQ/$PROJECT/*/Sample_${NAME}_*/*.fastq.gz -t 2 -o $DIR_FASTQC


# Trimming
DIR_TRIMMING=$DIR_MAP/$PROJECT/trimming
DIR_TRIMMING_FASTQC=$DIR_TRIMMING/fastqc
mkdir -p $DIR_TRIMMING_FASTQC
R1=$DIR_FASTQ/$PROJECT/*/Sample_${NAME}_*/*_R1_001.fastq.gz
R2=$DIR_FASTQ/$PROJECT/*/Sample_${NAME}_*/*_R2_001.fastq.gz
trim_galore --fastqc --fastqc_args "--outdir $DIR_TRIMMING_FASTQC" \
  --nextseq 20 --illumina --paired --length 15 $R1 $R2 -o $DIR_TRIMMING


# Mapping with bowtie2
mkdir -p $TMP/$NAME
cd $TMP/$NAME
R1=$DIR_TRIMMING/${NAME}_*_R1_001_val_1.fq.gz
R2=$DIR_TRIMMING/${NAME}_*_R2_001_val_2.fq.gz
bowtie2 -p $CPU -N 1 -x $GENOME -1 $R1 -2 $R2 -X 1000 |\
  samtools view -bS - > $NAME.bam


# Sort by coordinate
samtools sort $NAME.bam -o $NAME.st.bam -@ $CPU -m $(( 500 * $MEM / $CPU ))"M"
samtools index $NAME.st.bam


# Mapping statistics
picard CollectInsertSizeMetrics I=$NAME.st.bam \
                                O=$NAME.insert_size_metrics.1k.txt \
                                H=$NAME.insert_size_histogram.1k.pdf \
                                M=0.5 \
                                HISTOGRAM_WIDTH=1000 \
                                USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
samtools stats $NAME.st.bam > $NAME.st.txt


# Remove duplicates with picard
picard SortSam I=$NAME.st.bam O=$NAME.qst.bam SORT_ORDER=queryname \
               USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
picard MarkDuplicates I=$NAME.qst.bam \
                      O=$NAME.rmdup.bam \
                      METRICS_FILE=$NAME.marked_dup_metrics.txt \
                      REMOVE_DUPLICATES=true \
                      USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
                      # VALIDATION_STRINGENCY=SILENT


# Allele-specific allignment with SNPsplit if specified
if [[ $GENOME =~ .N-masked$ ]]; then
  SNP_FILE=$DIR_SNPSPLIT/$GENOME/all_*.txt.gz
 SNPsplit --snp_file $SNP_FILE --paired --conflicting \
   $NAME.rmdup.bam
  rm $NAME.rmdup.sortedByName.bam
  AS=( "" .genome1 .genome2 .unassigned .conflicting )
else
  AS=( "" )
fi


# Genome coordinate
GENOME_COORD=${GENOME%%.*}


# Get tag count maps
for I in "${!AS[@]}"
do
  # Extract Uniquely and properly mapped reads with samtools
  samtools view -u -f 2 -q 20 -h $NAME.rmdup${AS[$I]}.bam |\
    samtools sort -@ $CPU -m $(( 500 * $MEM / $CPU ))"M" \
      -o ${NAME}${AS[$I]}.uniq.bam
  samtools index ${NAME}${AS[$I]}.uniq.bam
 
 
  # Mapping statistics for uniquely mapped reads
  samtools idxstats ${NAME}${AS[$I]}.uniq.bam \
    > ${NAME}${AS[$I]}.uniq.idxstats.txt
  samtools stats ${NAME}${AS[$I]}.uniq.bam > ${NAME}${AS[$I]}.uniq.stats.txt
  picard CollectInsertSizeMetrics I=${NAME}${AS[$I]}.uniq.bam \
                            O=${NAME}${AS[$I]}.uniq.insert_size_metrics.1k.txt \
                            H=${NAME}${AS[$I]}.insert_size_histogram.1k.pdf \
                            M=0.5 \
                            HISTOGRAM_WIDTH=1000 \
                            USE_JDK_DEFLATER=true USE_JDK_INFLATER=true


  # Map resection endpoints
  echo "track type=wiggle_0" > ${NAME}${AS[$I]}.F.wig
  echo "track type=wiggle_0" > ${NAME}${AS[$I]}.R.wig

  if [ $ADAPTOR_SWAP = "FALSE" ]; then
    # Forward strand (99, 163)
    samtools view -f 99 ${NAME}${AS[$I]}.uniq.bam |\
     awk -v OFS="\t" '{print $3, $4}' |\
     sort -k1,1 -k2,2n |\
     uniq -c |\
     awk -v OFS="" \
       'BEGIN{CHR=""}
       {if(CHR!=$2){CHR=$2; print "variableStep chrom=", $2, " span=1"}}
       {print $3, "\t", $1}' \
       >> ${NAME}${AS[$I]}.F.wig

    # Reverse strand (83, 147)
    samtools view -f 83 ${NAME}${AS[$I]}.uniq.bam |\
     awk -v OFS="\t" '{print $3, $8 - $9 - 1}' |\
     sort -k1,1 -k2,2n |\
     uniq -c |\
     awk -v OFS="" \
       'BEGIN{CHR=""}
       {if(CHR!=$2){CHR=$2; print "variableStep chrom=", $2, " span=1"}}
       {print $3, "\t", $1}' \
       >> ${NAME}${AS[$I]}.R.wig
  elif [ $ADAPTOR_SWAP = "TRUE" ]; then
    # Forward strand (99, 163)
    samtools view -f 99 ${NAME}${AS[$I]}.uniq.bam |\
     awk -v OFS="\t" '{print $3, $4 + $9 - 1}' |\
     sort -k1,1 -k2,2n |\
     uniq -c |\
     awk -v OFS="" \
       'BEGIN{CHR=""}
       {if(CHR!=$2){CHR=$2; print "variableStep chrom=", $2, " span=1"}}
       {print $3, "\t", $1}' \
       >> ${NAME}${AS[$I]}.R.wig

    # Reverse strand (83, 147)
    samtools view -f 83 ${NAME}${AS[$I]}.uniq.bam |\
     awk -v OFS="\t" '{print $3, $8}' |\
     sort -k1,1 -k2,2n |\
     uniq -c |\
     awk -v OFS="" \
       'BEGIN{CHR=""}
       {if(CHR!=$2){CHR=$2; print "variableStep chrom=", $2, " span=1"}}
       {print $3, "\t", $1}' \
     >> ${NAME}${AS[$I]}.F.wig
  fi

  wigToBigWig ${NAME}${AS[$I]}.F.wig \
    $DIR_CHROM_SIZES/$GENOME_COORD.chrom.sizes ${NAME}${AS[$I]}.F.bw
  wigToBigWig ${NAME}${AS[$I]}.R.wig \
    $DIR_CHROM_SIZES/$GENOME_COORD.chrom.sizes ${NAME}${AS[$I]}.R.bw
  cat $NAME.F.wig $NAME.R.wig |\
    grep ^[1-9] |\
    awk -v OFS="\t" -v NAME=$NAME '{m+=$2} END {print NAME, m;}' \
    > ${NAME}${AS[$I]}.rpm.txt


  # RPM normalization and binning
  TOTAL_READS=`cat $NAME.rpm.txt | cut -f2`
  WIN=25
 
  for ST in F R
  do
    bigWigToBedGraph ${NAME}${AS[$I]}.$ST.bw ${NAME}${AS[$I]}.$ST.bedgraph
    cat ${NAME}${AS[$I]}.$ST.bedgraph |\
      awk -v OFS="\t" \
        '{for(i=$2+1; i<=$3; i++){for(j=1; j<=$4; j++){print $1, i-1, i}}}' \
      > ${NAME}${AS[$I]}.$ST.bed
    igvtools count -w $WIN ${NAME}${AS[$I]}.$ST.bed stdout \
      $DIR_CHROM_SIZES/$GENOME_COORD.chrom.sizes |\
      sed -n '/track type=wiggle_0/,$p' |\
      awk -v OFS="\t" -v TR=$TOTAL_READS \
        '{if($1 ~ /^[1-9]/){print $1, $2 * 1e6 / TR}else{print}}' \
      > ${NAME}${AS[$I]}.BIN.$WIN.bp.rpm.$ST.wig
    wigToBigWig ${NAME}${AS[$I]}.BIN.$WIN.bp.rpm.$ST.wig \
      $DIR_CHROM_SIZES/$GENOME_COORD.chrom.sizes \
      ${NAME}${AS[$I]}.BIN.$WIN.bp.rpm.$ST.bw
  done
done


# Get tag matrixes around peaks
WIN_L=5600
WIN_R=5601
DIR_PEAK=$DIR_PEAKS/$GENOME/primary
for PEAK in $DIR_PEAK/*.bed
do
  if [[ $PEAK =~ \*.bed$ ]]; then
    echo "No peak bed files in $DIR_PEAK"
    break
  fi
  for I in "${!AS[@]}"
  do
    for BW in ${NAME}${AS[$I]}.F.bw ${NAME}${AS[$I]}.R.bw
    do
      echo "-------------------------------------"
      echo "Date        : $( date )"
      echo "Peak        : ${PEAK##*/}"
      echo "Signal      : $BW"
      echo "Half window : $WIN_L"
      echo
      MAT=${PEAK##*/}.${BW%.bw}.halfwin.$WIN_L.matrix.txt.gz
      computeMatrix reference-point --referencePoint TSS -S $BW -R "$PEAK" \
        -b $WIN_L -a $WIN_R -out $MAT \
        --binSize 1 --sortRegions keep -p $CPU
      R --vanilla --slave --args $PEAK $MAT \
        < "$SDIR"/lib/matrix2Rdata.20191228.R
    done
    for BW in ${NAME}${AS[$I]}.F.bw
    do
      MAT_R_F=${PEAK##*/}.${BW%.bw}.halfwin.$WIN_L.matrix.Rdata
      DIR_MAT_R_CTRL=$DIR_MATS/control/$GENOME
      Rscript --vanilla --slave \
        "$SDIR"/lib/pdf.plot.mean.with.matrixRdata.20191228.R \
        $PEAK $MAT_R_F $DIR_MAT_R_CTRL
    done
  done
done


# Move the temporary directory to the MAP directory
DIR_BOWTIE2=$DIR_MAP/$PROJECT/$GENOME
mkdir -p $DIR_BOWTIE2
rm $NAME.bam $NAME.qst.bam $NAME.rmdup.bam $NAME.*.wig \
  $NAME.*.bed $NAME.*.bedgraph
mv $TMP/$NAME $DIR_BOWTIE2
rm $DIR_TRIMMING/*.fq.gz
