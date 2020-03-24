## ngspipe
pipeline for analysis of NGS data. From fastqs to BAM file


Usage:
./ngspipe.v1.sh FASTQ_R1 FASTQ_R2 SAMPLEID


Included steps:
- fastq clumping (to speed up and improve compression)
- fastqc and multiqc
- fastp (adapter trimming and quality filtering)
- bwa-mem mapping
- post aligment qc
- alignment around indels (GATK)
- base recalibrator (GATK BQRS)
- final qc on bam file