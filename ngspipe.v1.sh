#!/bin/bash

# create ngspipe conda environment
# tools installed: gatk bbmap samtools bwa openjdk fastqc multiqc picard libiconv r-gplots r-kernsmooth


conda activate ngspipe


FQ1=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R1_001.fastq.gz
FQ2=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R2_001.fastq.gz
FQ=ICP65
SUBSEQ=false

# set N threads
NT=16

# REEFERENCES
BWAGENOME=/mnt/jbod/nando/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
FASTAGENOME=/mnt/jbod/nando/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
KNOWN1=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.vcf.gz
KNOWN2=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
KNOWN3=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
FEATUREFILE=/home/ferdy/2020_ngs_design.bed

mkdir -p $FQ
cd $FQ

mkdir -p qc

# clumpify sort and speed up with compression
clumpify.sh in=$FQ1 in2=$FQ2 out=01.$FQ.R1.fq.gz out2=01.$FQ.R2.fq.gz reorder=p dedupe=f -Xmx8g  # set dedupe=t if willing to deduplicate, but why?


mkdir -p qc/01qc
fastqc -t $NT -o qc/01qc 01.$FQ.R1.fq.gz 01.$FQ.R2.fq.gz
multiqc -f -n 00.pre.trimming -o qc/ qc/01qc/


# remove low quality reads
if [ SUBSEQ = true ]; then
  echo "Subseq is active"
  # subseq
  seqtk sample -s100 01.$FQ.R1.fq.gz 0.1 | pigz -p $NT > 01.$FQ.R1.sub.fq.gz
  seqtk sample -s100 01.$FQ.R2.fq.gz 0.1 | pigz -p $NT > 01.$FQ.R2.sub.fq.gz
  
  fastp \
--in1 01.$FQ.R1.sub.fq.gz \
--in2 01.$FQ.R2.sub.fq.gz \
--cut_front \
--cut_tail \
--thread $NT \
--detect_adapter_for_pe \
--out1 02.$FQ.R1.fq.gz \
--out2 02.$FQ.R2.fq.gz \
--json qc/trim.$FQ.fastp.json \
--html qc/trim.$FQ.fastp.html > qc/$FQ.fastp.out 2>&1

  else
  echo "Subseq not active"

  fastp \
--in1 01.$FQ.R1.fq.gz \
--in2 01.$FQ.R2.fq.gz \
--cut_front \
--cut_tail \
--thread $NT \
--detect_adapter_for_pe \
--out1 02.$FQ.R1.fq.gz \
--out2 02.$FQ.R2.fq.gz \
--json qc/trim.$FQ.fastp.json \
--html qc/trim.$FQ.fastp.html > qc/$FQ.fastp.out 2>&1

fi


multiqc -m fastp -f -n 01.post.trimming.fastp -o qc/ qc/*fastp.json
rm 01.$FQ* qc/trim.$FQ.fastp.*



# align with BWA MEM
bwa mem -t $NT -R @RG\\tID:$FQ\\tSM:$FQ\\tPL:ILLUMINA\\tLB:$FQ $BWAGENOME 02.$FQ.R1.fq.gz 02.$FQ.R2.fq.gz | samtools sort -@ $NT -O BAM -o 03.$FQ.bam -
rm 01.$FQ* 02.$FQ*

# mark DUPS
picard MarkDuplicates \
I=03.$FQ.bam \
O=04.$FQ.marked.bam \
M=qc/$FQ.sorted.picard.metrics \
REMOVE_DUPLICATES=false \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=./

rm 03.$FQ.bam

multiqc -m picard -f -n 02.picard.markdups -o qc/ qc/$FQ.sorted.picard.metrics
rm qc/$FQ.sorted.picard.metrics

# realign around indels
parallel -j $NT "
# generate chromosome interval
gatk3 -T RealignerTargetCreator -R $FASTAGENOME -I 04.$FQ.marked.bam -o $FQ.chr{}.intervals -L chr{}
# local realignment
gatk3 -T IndelRealigner -R $FASTAGENOME -I 04.$FQ.marked.bam -targetIntervals $FQ.chr{}.intervals -L chr{} -o $FQ.chr{}.realigned.bam
" ::: {1..22} X Y M

# concatenate realigned reads and unmapped ones
samtools cat $FQ.chr*.realigned.bam | samtools sort -@$NT > 05.$FQ.realigned.bam
samtools index -@$NT 05.$FQ.realigned.bam
rm $FQ.chr*.realigned.ba* $FQ.chr*.intervals 04.$FQ* $FQ.recal_data.table


# recalibrate reads NB: add other -knownSites for additional vcfs

gatk3 \
-nct $NT \
-T BaseRecalibrator \
-R $FASTAGENOME \
-I 05.$FQ.realigned.bam \
-knownSites $KNOWN1 \
-knownSites $KNOWN2 \
-knownSites $KNOWN3 \
-o $FQ.recal_data.table


# apply recalibration
gatk3 \
-nct $NT \
-T PrintReads \
-R $FASTAGENOME \
-I 05.$FQ.realigned.bam \
-BQSR $FQ.recal_data.table \
-o $FQ.final.bam


rm $FQ.recal_data.table 05.$FQ.realigned.bam


# quality control final bam

qualimap bamqc \
-bam $FQ.final.bam \
--feature-file $FEATUREFILE \
-nt $NT \
-sd \
-outdir qc/ \
-outformat PDF:HTML

multiqc -m qualimap -f -n 03.bamqc.report.final -o qc/ qc/raw_data_qualimapReport/ qc/genome_results.txt
rm -r qc/$FQ.sorted.picard.metrics qc/trim.$FQ*

#######################
## END
#######################