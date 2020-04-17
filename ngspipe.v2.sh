#!/bin/bash


# 20200324 - meeting. Step 1 validation outputs: coverage and variant called with old and new bam files

# # if miniconda or conda is not installed, install it with:
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# # create a conda environment with all the needed tools (NB gatk version 3)
# conda create -n ngspipe gatk bbmap samtools bwa openjdk fastqc multiqc picard fastuniq libiconv r-gplots r-kernsmooth qualimap fastp seqtk parallel -y
# # register gatk version
# conda activate ngspipe
# gatk3-register /mnt/jbod/common/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
# # remove parallel citation message
# parallel --citation
# conda deactivate
# # now the environment is ready to go

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
# FQ=ICP65.GRCh37

###########
## START ##
###########
echo -e "\e[33m[NGSPIPE] Activating conda environment \e[39m"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ngspipe
SUBSEQ="false"
NT=8 # set N threads

# REEFERENCES
# please use always GRCh37 and download with aws from igenomes (https://ewels.github.io/AWS-iGenomes/)
# aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/ /mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/
BWAGENOME=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Sequence/BWAIndex/human_g1k_v37_decoy.fasta
FASTAGENOME=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta
KNOWN1=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf
KNOWN2=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf
KNOWN3=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf 
FEATUREFILE=/mnt/jbod/common/Homo_sapiens/2020_ngs_design.bed

mkdir -p $FQ
cd $FQ

mkdir -p qc

# clumpify sort and speed up with compression
echo -e "\e[33m[NGSPIPE] Clumping reads \e[39m"
clumpify.sh in=$FQ1 in2=$FQ2 out=01.$FQ.R1.fq.gz out2=01.$FQ.R2.fq.gz reorder=p dedupe=f -Xmx8g > /dev/null 2>&1

echo -e "\e[33m[NGSPIPE] raw reads QC \e[39m"
mkdir -p qc/01qc
fastqc -t $NT -o qc/01qc 01.$FQ.R1.fq.gz 01.$FQ.R2.fq.gz 
multiqc -f -n 00.pre.trimming -o qc/ qc/01qc/ > /dev/null 2>&1

# remove low quality reads
echo -e "\e[33m[NGSPIPE] Trimming and filtering reads \e[39m"
if [ "$SUBSEQ" = true ]; then
  echo "[NGSPIPE] Subseq is active, only 10% of reads will be used"
  # subseq
  echo "[NGSPIPE] Subsampling reads"
  seqtk sample -s100 01.$FQ.R1.fq.gz 0.1 | pigz -p $NT > 01.$FQ.R1.sub.fq.gz
  seqtk sample -s100 01.$FQ.R2.fq.gz 0.1 | pigz -p $NT > 01.$FQ.R2.sub.fq.gz
  echo "[NGSPIPE] Trimming and filtering reads"
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
  echo "[NGSPIPE] Subseq is not active, 100% of reads will be used"
  echo "[NGSPIPE] Trimming and filtering reads"
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

echo -e "\e[33m[NGSPIPE] Post-filtering QC \e[39m"
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
# -maxReads default is 20000
# add known indels to the model
echo -e "\e[33m[NGSPIPE] Realign around indels \e[39m"
parallel -j $NT "
	# generate chromosome interval
	gatk3 -T RealignerTargetCreator -R $FASTAGENOME -I 04.$FQ.marked.bam -o $FQ.chr{}.intervals -L {}
	# local realignment
	gatk3 -T IndelRealigner -R $FASTAGENOME -I 04.$FQ.marked.bam -maxReads 50000 -known $KNOWN2 -targetIntervals $FQ.chr{}.intervals -L {} -o $FQ.chr{}.realigned.bam
" ::: {1..22} X Y MT

# concatenate realigned reads and unmapped ones
samtools cat $FQ.chr*.realigned.bam | samtools sort -@$NT > 05.$FQ.realigned.bam
samtools index -@$NT 05.$FQ.realigned.bam
rm $FQ.chr*.realigned.ba* $FQ.chr*.intervals 04.$FQ* $FQ.recal_data.table


# recalibrate reads 
echo -e "\e[33m[NGSPIPE] Calculating recalibration scores (GATK BQSR) \e[39m"
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
echo -e "\e[33m[NGSPIPE] Recalibrating (GATK BQSR) \e[39m"
gatk3 \
-nct $NT \
-T PrintReads \
-R $FASTAGENOME \
-I 05.$FQ.realigned.bam \
-BQSR $FQ.recal_data.table \
-o $FQ.prefinal.bam


rm $FQ.recal_data.table 05.$FQ.realigned.bam 05.$FQ.realigned.bam.bai
samtools sort -@$NT $FQ.prefinal.bam > $FQ.final.bam
samtools index $FQ.final.bam
rm $FQ.prefinal.ba*

# fix feature file
sed 's/chrM/MT/g' $FEATUREFILE > featurefile.bed
sed -i 's/chr//g' featurefile.bed

# quality control final bam
echo -e "\e[33m[NGSPIPE] Quality control on final bam \e[39m"
qualimap bamqc \
-bam $FQ.final.bam \
--feature-file featurefile.bed \
-nt $NT \
-sd \
-outdir qc/ \
-outformat PDF:HTML

multiqc -m qualimap -f -n 03.bamqc.report.final -o qc/ qc/raw_data_qualimapReport/ qc/genome_results.txt
rm featurefile.bed

cd ..
conda deactivate

#########
## END ##
#########
