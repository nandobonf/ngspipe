# ngspipe
pipeline for analysis of NGS data. From fastqs to BAM file


## Usage

```bash
./ngspipe.v1.sh {FASTQ_R1} {FASTQ_R2} {SAMPLEID}
```

## Included steps:
- fastq clumping (to speed up and improve compression)
- fastqc and multiqc (raw fastq plots and visualization)
- fastp (adapter trimming and quality filtering)
- bwa-mem mapping to GRCh37 assembly
- post aligment qc
- alignment around indels (GATK)
- base recalibrator (GATK BQSR)
- final qc on bam file

## In order to install the conda environment 
```bash
# if miniconda or conda is not installed, install it with:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# create a conda environment with all the needed tools (NB gatk version 3)
conda create -n ngspipe gatk bbmap samtools bwa openjdk fastqc multiqc picard fastuniq libiconv r-gplots r-kernsmooth qualimap fastp seqtk parallel -y
# register gatk version, download from gatk website the file "GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2"
conda activate ngspipe
gatk3-register /mnt/jbod/common/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
# remove parallel citation message
parallel --citation
conda deactivate
```

## before running the pipeline, download the required files and change their variables accordingly
```bash
# REEFERENCES
# please use always GRCh37 (GATK bundle) and download with aws from igenomes (https://ewels.github.io/AWS-iGenomes/)
# aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/ /mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/
# path the the BWA genome index
BWAGENOME=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Sequence/BWAIndex/human_g1k_v37_decoy.fasta
# path the the GRCh37 fasta file
FASTAGENOME=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta
# list of VCF files ok known variants for base recalibration
KNOWN1=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf
KNOWN2=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf
KNOWN3=/mnt/jbod/nando/GATK.bundle/references/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf 
# bed file of target regions for targeted panel sequencing
FEATUREFILE=/mnt/jbod/common/Homo_sapiens/2020_ngs_design.bed

```
