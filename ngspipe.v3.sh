#!/bin/bash

# 20210304 implemente gatk4

# # if miniconda or conda is not installed, install it with:
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# # create a conda environment with all the needed tools (NB gatk version 3)
# conda create -n ngspipe ascli gatk4 bbmap samtools bwa openjdk fastqc multiqc picard fastuniq libiconv r-gplots r-kernsmooth qualimap fastp seqtk parallel sambamba -y
# conda activate ngspipe
# # remove parallel citation message
# parallel --citation
# conda deactivate
# # now the environment is ready to go

# usage:
#~/bin/ngspipeline_v01_FB.sh /mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R1_001.fastq.gz /mnt/jbod/common/xNando/Run_Wei/ICP65_S14_L001_R2_001.fastq.gz ICP65

FQ1=$1
FQ2=$2
FQ=$3

# # # ONLY FOR TESTING
# FQ1=/srv/ngsdata/RAWDATA/WES_V7_germline_NB_biodiversa/1740-74/1740D-74-09_S62_L003_R1_001.fastq.gz
# FQ2=/srv/ngsdata/RAWDATA/WES_V7_germline_NB_biodiversa/1740-74/1740D-74-09_S62_L003_R2_001.fastq.gz
# FQ=MC

# FQ1=/srv/ngsdata/nando/test.R1.fq.gz
# FQ2=/srv/ngsdata/nando/test.R2.fq.gz
# FQ=TEST02


echo -e "\e[33m[NGSPIPE] FASTQ1 is $FQ1 \e[39m"
echo -e "\e[33m[NGSPIPE] FASTQ2 is $FQ2 \e[39m"
echo -e "\e[33m[NGSPIPE] SAMPLE NAME is $FQ \e[39m"


###########
## START ##
###########
echo -e "\e[33m[NGSPIPE] Activating conda environment \e[39m"
source /srv/ngsdata/nando/anaconda3/etc/profile.d/conda.sh
conda activate ngspipev2
NT=8 # set N threads
export TEMPDIR=/srv/ngsdata/nando/tmp
export _JAVA_OPTIONS=-Djava.io.tmpdir=/srv/ngsdata/nando/tmp


# REEFERENCES
# please use always GRCh37 and download with aws from igenomes (https://ewels.github.io/AWS-iGenomes/)
# aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/ /mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/
BWAGENOME=/srv/ngsdata/nando/references/Homo_sapiens/GATK/GRCh37/Sequence/BWAIndex/human_g1k_v37_decoy.fasta
FASTAGENOME=/srv/ngsdata/nando/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta
KNOWN1=/srv/ngsdata/nando/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf
KNOWN2=/srv/ngsdata/nando/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf
KNOWN3=/srv/ngsdata/nando/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf 
#FEATUREFILE=/mnt/jbod/common/Homo_sapiens/2020_ngs_design.bed

mkdir -p $FQ
cd $FQ

mkdir -p qc


# remove low quality reads
echo -e "\e[33m[NGSPIPE] Trimming and filtering reads \e[39m"

fastp \
--in1 $FQ1 \
--in2 $FQ2 \
--cut_front \
--cut_tail \
--thread $NT \
--detect_adapter_for_pe \
--out1 02.$FQ.R1.fq.gz \
--out2 02.$FQ.R2.fq.gz \
--json qc/trim.$FQ.fastp.json \
--html qc/trim.$FQ.fastp.html > qc/$FQ.fastp.out 2>&1


echo -e "\e[33m[NGSPIPE] Post-filtering QC \e[39m"
multiqc -m fastp -f -n 01.post.trimming.fastp -o qc/ qc/*fastp.json
rm qc/trim.$FQ.fastp.*


# align with BWA MEM
echo -e "\e[33m[NGSPIPE] Aligning reads \e[39m"
bwa mem -Ma -t $NT -R @RG\\tID:$FQ\\tSM:$FQ\\tPL:ILLUMINA\\tLB:$FQ $BWAGENOME 02.$FQ.R1.fq.gz 02.$FQ.R2.fq.gz | samtools sort -@ $NT -O BAM -o 03.$FQ.bam -
rm 02.$FQ*


echo -e "\e[33m[NGSPIPE] Marking duplicates (MarkDuplicates)\e[39m"
gatk --java-options "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=8" \
MarkDuplicates \
-I 03.$FQ.bam \
-O 04.$FQ.marked.bam \
-M qc/$FQ.sorted.picard.metrics \
--VALIDATION_STRINGENCY LENIENT \
--CREATE_INDEX TRUE \
--TMP_DIR /srv/ngsdata/nando/tmp



rm 03.$FQ.bam

multiqc -m picard -f -n 02.picard.markdups -o qc/ qc/$FQ.sorted.picard.metrics
rm qc/$FQ.sorted.picard.metrics



echo -e "\e[33m[NGSPIPE] Calculating recalibration scores (GATK BaseRecalibrator) \e[39m"
gatk --java-options "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=8" \
BaseRecalibrator \
-R $FASTAGENOME \
-I 04.$FQ.marked.bam \
-known-sites $KNOWN1 \
-known-sites $KNOWN2 \
-known-sites $KNOWN3 \
-O $FQ.recal_data.table



# apply recalibration (Spark)
echo -e "\e[33m[NGSPIPE] Recalibrating (GATK ApplyBQSRS) \e[39m"
gatk --java-options "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=8" \
ApplyBQSR \
-R $FASTAGENOME \
-I 04.$FQ.marked.bam \
-bqsr $FQ.recal_data.table \
--static-quantized-quals 10 \
--static-quantized-quals 20 \
--static-quantized-quals 30 \
-O $FQ.bam

mv $FQ.bai $FQ.bam.bai

# quality control final bam
echo -e "\e[33m[NGSPIPE] Quality control on final bam \e[39m"
qualimap bamqc \
-bam $FQ.bam \
-nt $NT \
-sd \
-outdir qc/ \
--java-mem-size=16G \
-outformat PDF:HTML

multiqc -m qualimap -f -n 03.bamqc.report.final -o qc/ qc/raw_data_qualimapReport/ qc/genome_results.txt
rm 04.$FQ.marked.b* $FQ.recal_data.table


cd ..
conda deactivate

echo -e "\e[33m[NGSPIPE] Analysis complete \e[39m"

#########
## END ##
#########

