## ngspipe
pipeline for analysis of NGS data. From fastqs to BAM file


## Usage

```bash
./ngspipe.v1.sh {FASTQ_R1} {FASTQ_R2} {SAMPLEID}
```

##Included steps:
- fastq clumping (to speed up and improve compression)
- fastqc and multiqc
- fastp (adapter trimming and quality filtering)
- bwa-mem mapping
- post aligment qc
- alignment around indels (GATK)
- base recalibrator (GATK BQRS)
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
conda create -n ngspipe gatk bbmap samtools bwa openjdk fastqc multiqc picard libiconv r-gplots r-kernsmooth qualimap fastp seqtk parallel -y
# register gatk version
conda activate ngspipe
gatk3-register /mnt/jbod/common/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
# remove parallel citation message
parallel --citation
conda deactivate
```