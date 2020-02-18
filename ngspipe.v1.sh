#!/bin/bash

# create ngspipe conda environment
# tools installed: gatk bbmap samtools bwa openjdk fastqc multiqc picard libiconv r-gplots r-kernsmooth
# conda

conda activate ngspipe

# set input files
$1=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R1_001.fastq.gz # R1
$2=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R2_001.fastq.gz # R2
$3=ICP65                                                          # sample name
$4=16                                                             # threads / cores

SUBSEQ=false


# REEFERENCES
BWAGENOME=/mnt/jbod/nando/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
FASTAGENOME=/mnt/jbod/nando/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
KNOWN1=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.vcf.gz
KNOWN2=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
KNOWN3=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
FEATUREFILE=/home/ferdy/2020_ngs_design.bed

mkdir -p $3
cd $3

mkdir -p qc

# clumpify sort and speed up with compression
clumpify.sh in=$1 in2=$2 out=01.$3.R1.fq.gz out2=01.$3.R2.fq.gz reorder=p dedupe=f -Xmx8g  # set dedupe=t if willing to deduplicate, but why?


mkdir -p qc/01qc
fastqc -t $4 -o qc/01qc 01.$3.R1.fq.gz 01.$3.R2.fq.gz
multiqc -f -n 00.pre.trimming -o qc/ qc/01qc/


# remove low quality reads
if [ SUBSEQ = true ]; then
  echo "Subseq is active"
  # subseq
  seqtk sample -s100 01.$3.R1.fq.gz 0.1 | pigz -p $4 > 01.$3.R1.sub.fq.gz
  seqtk sample -s100 01.$3.R2.fq.gz 0.1 | pigz -p $4 > 01.$3.R2.sub.fq.gz
  
  fastp \
--in1 01.$3.R1.sub.fq.gz \
--in2 01.$3.R2.sub.fq.gz \
--cut_front \
--cut_tail \
--thread $4 \
--detect_adapter_for_pe \
--out1 02.$3.R1.fq.gz \
--out2 02.$3.R2.fq.gz \
--json qc/trim.$3.fastp.json \
--html qc/trim.$3.fastp.html > qc/$3.fastp.out 2>&1

  else
  echo "Subseq not active"

  fastp \
--in1 01.$3.R1.fq.gz \
--in2 01.$3.R2.fq.gz \
--cut_front \
--cut_tail \
--thread $4 \
--detect_adapter_for_pe \
--out1 02.$3.R1.fq.gz \
--out2 02.$3.R2.fq.gz \
--json qc/trim.$3.fastp.json \
--html qc/trim.$3.fastp.html > qc/$3.fastp.out 2>&1

fi


multiqc -m fastp -f -n 01.post.trimming.fastp -o qc/ qc/*fastp.json
rm 01.$3* qc/trim.$3.fastp.*



# align with BWA MEM
bwa mem -t $4 -R @RG\\tID:$3\\tSM:$3\\tPL:ILLUMINA\\tLB:$3 $BWAGENOME 02.$3.R1.fq.gz 02.$3.R2.fq.gz | samtools sort -@ $4 -O BAM -o 03.$3.bam -
rm 01.$3* 02.$3*

# mark DUPS
picard MarkDuplicates \
I=03.$3.bam \
O=04.$3.marked.bam \
M=qc/$3.sorted.picard.metrics \
REMOVE_DUPLICATES=false \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=./

rm 03.$3.bam

multiqc -m picard -f -n 02.picard.markdups -o qc/ qc/$3.sorted.picard.metrics
rm qc/$3.sorted.picard.metrics

# realign around indels
parallel -j $4 "
# generate chromosome interval
gatk3 -T RealignerTargetCreator -R $FASTAGENOME -I 04.$3.marked.bam -o $3.chr{}.intervals -L chr{}
# local realignment
gatk3 -T IndelRealigner -R $FASTAGENOME -I 04.$3.marked.bam -targetIntervals $3.chr{}.intervals -L chr{} -o $3.chr{}.realigned.bam
" ::: {1..22} X Y M

# concatenate realigned reads and unmapped ones
samtools cat $3.chr*.realigned.bam | samtools sort -@$4 > 05.$3.realigned.bam
samtools index -@$4 05.$3.realigned.bam
rm $3.chr*.realigned.ba* $3.chr*.intervals 04.$3* $3.recal_data.table


# recalibrate reads NB: add other -knownSites for additional vcfs

gatk3 \
-nct $4 \
-T BaseRecalibrator \
-R $FASTAGENOME \
-I 05.$3.realigned.bam \
-knownSites $KNOWN1 \
-knownSites $KNOWN2 \
-knownSites $KNOWN3 \
-o $3.recal_data.table


# apply recalibration
gatk3 \
-nct $4 \
-T PrintReads \
-R $FASTAGENOME \
-I 05.$3.realigned.bam \
-BQSR $3.recal_data.table \
-o $3.final.bam


rm $3.recal_data.table 05.$3.realigned.bam


# quality control final bam

qualimap bamqc \
-bam $3.final.bam \
--feature-file $FEATUREFILE \
-nt $4 \
-sd \
-outdir qc/ \
-outformat PDF:HTML

multiqc -m qualimap -f -n 03.bamqc.report.final -o qc/ qc/raw_data_qualimapReport/ qc/genome_results.txt
rm -r qc/$3.sorted.picard.metrics qc/trim.$3*

#######################
## END
#######################