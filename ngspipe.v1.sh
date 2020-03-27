#!/bin/bash


# 20200324 - meeting. Step 1 validation outputs: coverage and variant called with old and new bam files

# # if miniconda or conda is not installed, install it with:
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# # create a conda environment with all the needed tools (NB gatk version 3)
# conda create -n ngspipe gatk bbmap samtools bwa openjdk fastqc multiqc picard libiconv r-gplots r-kernsmooth qualimap fastp seqtk parallel -y
# # register gatk version
# conda activate ngspipe
# gatk3-register /mnt/jbod/common/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
# # remove parallel citation message
# parallel --citation
# conda deactivate

# usage:
#~/bin/ngspipeline_v01_FB.sh /mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R1_001.fastq.gz /mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R2_001.fastq.gz ICP65

FQ1=$1
FQ2=$2
FQ=$3

echo -e "\e[33m[NGSPIPE] FASTQ1 is $FQ1 \e[39m"
echo -e "\e[33m[NGSPIPE] FASTQ2 is $FQ2 \e[39m"
echo -e "\e[33m[NGSPIPE] SAMPLE NAME is $FQ \e[39m"

# # ONLY FOR TESTING
# FQ1=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R1_001.fastq.gz
# FQ2=/mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R2_001.fastq.gz
# FQ=ICP65

###########
## START ##
###########
echo -e "\e[33m[NGSPIPE] Activating conda environment \e[39m"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ngspipe
SUBSEQ=false
NT=8 # set N threads

# REEFERENCES
BWAGENOME=/mnt/jbod/common/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
FASTAGENOME=/mnt/jbod/common/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
KNOWN1=/mnt/jbod/common/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.vcf.gz
KNOWN2=/mnt/jbod/common/Homo_sapiens/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
KNOWN3=/mnt/jbod/common/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
FEATUREFILE=/mnt/jbod/common/Homo_sapiens/2020_ngs_design.bed

mkdir -p $FQ
cd $FQ

mkdir -p qc

# clumpify sort and speed up with compression
echo -e "\e[33m[NGSPIPE] Clumping reads \e[39m"
clumpify.sh in=$FQ1 in2=$FQ2 out=01.$FQ.R1.fq.gz out2=01.$FQ.R2.fq.gz reorder=p dedupe=f -Xmx8g 

mkdir -p qc/01qc
fastqc -t $NT -o qc/01qc 01.$FQ.R1.fq.gz 01.$FQ.R2.fq.gz
multiqc -f -n 00.pre.trimming -o qc/ qc/01qc/

# remove low quality reads
echo -e "\e[33m[NGSPIPE] Trimming and filtering reads \e[39m"
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
echo -e "\e[33m[NGSPIPE] Aligning reads \e[39m"
bwa mem -t $NT -R @RG\\tID:$FQ\\tSM:$FQ\\tPL:ILLUMINA\\tLB:$FQ $BWAGENOME 02.$FQ.R1.fq.gz 02.$FQ.R2.fq.gz | samtools sort -@ $NT -O BAM -o 03.$FQ.bam -
rm 02.$FQ*

# mark DUPS
echo -e "\e[33m[NGSPIPE] Marking duplicates \e[39m"
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
echo -e "\e[33m[NGSPIPE] Realign around indels \e[39m"
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


# recalibrate reads 
echo -e "\e[33m[NGSPIPE] Recalibrate alignment scores (GATK BQSR) \e[39m"
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


rm $FQ.recal_data.table 05.$FQ.realigned.bam 05.$FQ.realigned.bam.bai

# quality control final bam
echo -e "\e[33m[NGSPIPE] Quality control on final bam \e[39m"
qualimap bamqc \
-bam $FQ.final.bam \
--feature-file $FEATUREFILE \
-nt $NT \
-sd \
-outdir qc/ \
-outformat PDF:HTML

multiqc -m qualimap -f -n 03.bamqc.report.final -o qc/ qc/raw_data_qualimapReport/ qc/genome_results.txt
rm -r qc/$FQ.sorted.picard.metrics qc/trim.$FQ*

cd ..
conda deactivate

#########
## END ##
#########
