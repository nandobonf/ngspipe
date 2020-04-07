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
- bwa-mem mapping to Hg19 assembly
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

## before running the pipeline, download the required files and change their link accordingly
```bash
# path the the BWA genome index
BWAGENOME=/mnt/jbod/common/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
# path the the Hg19 fasta file
FASTAGENOME=/mnt/jbod/common/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
# list of VCF files ok known variants for base recalibration
KNOWN1=/mnt/jbod/common/Homo_sapiens/GATK/hg19/dbsnp_138.hg19.vcf.gz
KNOWN2=/mnt/jbod/common/Homo_sapiens/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
KNOWN3=/mnt/jbod/common/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
# bed file of the target regions
FEATUREFILE=/mnt/jbod/common/Homo_sapiens/2020_ngs_design.bed
```
