#!/bin/bash

# create ngspipe conda environment
# tools installed: gatk bbmap samtools bwa openjdk fastqc multiqc picard libiconv r-gplots r-kernsmooth
# conda

conda activate ngspipe

# set input files
FQ1=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R1_001.fastq.gz # R1
FQ2=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R2_001.fastq.gz # R2
NAME=ICP65                                                         # sample name
NT=16                                                              # threads / cores
SUBSEQ=false                                                       # subsample 10% of the reads (only for testing)


# REEFERENCES
BWAGENOME=/mnt/jbod/nando/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
FASTAGENOME=/mnt/jbod/nando/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
KNOWN1=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.vcf.gz
KNOWN2=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
KNOWN3=/mnt/jbod/nando/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
FEATUREFILE=/home/ferdy/2020_ngs_design.bed

mkdir -p $NAME
cd $NAME

mkdir -p qc

# clumpify sort and speed up with compression
clumpify.sh in=$FQ1 in2=$FQ2 out=01.$NAME.R1.fq.gz out2=01.$NAME.R2.fq.gz reorder=p dedupe=f -Xmx8g  # set dedupe=t if willing to deduplicate, but why?

mkdir -p qc/01qc
fastqc -t $NT -o qc/01qc 01.$NAME.R1.fq.gz 01.$NAME.R2.fq.gz
multiqc -f -n 00.pre.trimming -o qc/ qc/01qc/

# remove low quality reads
if [ SUBSEQ = true ]; then
  echo "Subseq is active"
  # subseq
  seqtk sample -s100 01.$NAME.R1.fq.gz 0.1 | pigz -p $NT > 01.$NAME.R1.sub.fq.gz
  seqtk sample -s100 01.$NAME.R2.fq.gz 0.1 | pigz -p $NT > 01.$NAME.R2.sub.fq.gz
  
  fastp \
  --in1 01.$NAME.R1.sub.fq.gz \
  --in2 01.$NAME.R2.sub.fq.gz \
  --cut_front \
  --cut_tail \
  --thread $NT \
  --detect_adapter_for_pe \
  --out1 02.$NAME.R1.fq.gz \
  --out2 02.$NAME.R2.fq.gz \
  --json qc/trim.$NAME.fastp.json \
  --html qc/trim.$NAME.fastp.html > qc/$NAME.fastp.out 2>&1

  else
  echo "Subseq not active"

  fastp \
  --in1 01.$NAME.R1.fq.gz \
  --in2 01.$NAME.R2.fq.gz \
  --cut_front \
  --cut_tail \
  --thread $NT \
  --detect_adapter_for_pe \
  --out1 02.$NAME.R1.fq.gz \
  --out2 02.$NAME.R2.fq.gz \
  --json qc/trim.$NAME.fastp.json \
  --html qc/trim.$NAME.fastp.html > qc/$NAME.fastp.out 2>&1

fi


multiqc -m fastp -f -n 01.post.trimming.fastp -o qc/ qc/*fastp.json
rm 01.$NAME* qc/trim.$NAME.fastp.*



# align with BWA MEM
bwa mem -t $NT -R @RG\\tID:$NAME\\tSM:$NAME\\tPL:ILLUMINA\\tLB:$NAME $BWAGENOME 02.$NAME.R1.fq.gz 02.$NAME.R2.fq.gz | samtools sort -@ $NT -O BAM -o 03.$NAME.bam -
rm 01.$NAME* 02.$NAME*

# mark DUPS
picard MarkDuplicates \
I=03.$NAME.bam \
O=04.$NAME.marked.bam \
M=qc/$NAME.sorted.picard.metrics \
REMOVE_DUPLICATES=false \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=./

rm 03.$NAME.bam

multiqc -m picard -f -n 02.picard.markdups -o qc/ qc/$NAME.sorted.picard.metrics
rm qc/$NAME.sorted.picard.metrics

# realign around indels
parallel -j $NT "
# generate chromosome interval
gatk3 -T RealignerTargetCreator -R $FASTAGENOME -I 04.$NAME.marked.bam -o $NAME.chr{}.intervals -L chr{}
# local realignment
gatk3 -T IndelRealigner -R $FASTAGENOME -I 04.$NAME.marked.bam -targetIntervals $NAME.chr{}.intervals -L chr{} -o $NAME.chr{}.realigned.bam
" ::: {1..22} X Y M

# concatenate realigned reads and unmapped ones
samtools cat $NAME.chr*.realigned.bam | samtools sort -@$NT > 05.$NAME.realigned.bam
samtools index -@$NT 05.$NAME.realigned.bam
rm $NAME.chr*.realigned.ba* $NAME.chr*.intervals 04.$NAME* $NAME.recal_data.table


# recalibrate reads NB: add other -knownSites for additional vcfs

gatk3 \
-nct $NT \
-T BaseRecalibrator \
-R $FASTAGENOME \
-I 05.$NAME.realigned.bam \
-knownSites $KNOWN1 \
-knownSites $KNOWN2 \
-knownSites $KNOWN3 \
-o $NAME.recal_data.table


# apply recalibration
gatk3 \
-nct $NT \
-T PrintReads \
-R $FASTAGENOME \
-I 05.$NAME.realigned.bam \
-BQSR $NAME.recal_data.table \
-o $NAME.final.bam


rm $NAME.recal_data.table 05.$NAME.realigned.bam


# quality control final bam

qualimap bamqc \
-bam $NAME.final.bam \
--feature-file $FEATUREFILE \
-nt $NT \
-sd \
-outdir qc/ \
-outformat PDF:HTML

multiqc -m qualimap -f -n 03.bamqc.report.final -o qc/ qc/raw_data_qualimapReport/ qc/genome_results.txt
rm -r qc/$NAME.sorted.picard.metrics qc/trim.$NAME*

#######################
## END
#######################